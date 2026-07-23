#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
})

parse_args <- function(args) {
  opts <- list(
    delta_file = "SNPGenie/analysis/delta_selection/delta_per_sample.tsv",
    outdir = "SNPGenie/analysis/delta_selection",
    manifest = "",
    write_detailed = FALSE
  )

  if (length(args) == 0) {
    return(opts)
  }

  for (arg in args) {
    if (grepl("^--delta-file=", arg)) {
      opts$delta_file <- sub("^--delta-file=", "", arg)
    } else if (grepl("^--outdir=", arg)) {
      opts$outdir <- sub("^--outdir=", "", arg)
    } else if (grepl("^--manifest=", arg)) {
      opts$manifest <- sub("^--manifest=", "", arg)
    } else if (arg == "--write-detailed") {
      opts$write_detailed <- TRUE
    }
  }
  opts
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

run_limma_for_threshold <- function(df_thr, threshold, outdir, manifest_df, write_detailed = FALSE) {
  products <- sort(unique(df_thr$product))
  samples <- sort(unique(df_thr$sample))

  matrix_delta <- matrix(NA_real_, nrow = length(products), ncol = length(samples))
  rownames(matrix_delta) <- products
  colnames(matrix_delta) <- samples

  for (i in seq_len(nrow(df_thr))) {
    rr <- df_thr[i, ]
    matrix_delta[rr$product, rr$sample] <- rr$delta_piN_minus_piS
  }

  keep_genes <- apply(matrix_delta, 1, function(x) all(!is.na(x)))
  matrix_delta <- matrix_delta[keep_genes, , drop = FALSE]

  # Map sample name to DPI using manifest
  get_dpi <- function(sample_name) {
    match_row <- manifest_df[manifest_df$bam_name == sample_name, ]
    if (nrow(match_row) > 0) {
      return(paste0("dpi", match_row$dpi[1]))
    }
    return(NA_character_)
  }

  sample_info <- data.frame(
    sample = colnames(matrix_delta),
    dpi = vapply(colnames(matrix_delta), get_dpi, character(1)),
    stringsAsFactors = FALSE
  )

  # Filter out samples with missing DPI
  sample_info <- sample_info[!is.na(sample_info$dpi), , drop = FALSE]
  
  if (nrow(sample_info) == 0) {
    message("No samples with valid DPI values found for threshold ", threshold)
    return(NULL)
  }

  unique_dpis <- sort(unique(sample_info$dpi))
  sample_info$dpi <- factor(sample_info$dpi, levels = unique_dpis)

  matrix_delta <- matrix_delta[, sample_info$sample, drop = FALSE]

  design <- model.matrix(~ 0 + dpi, data = sample_info)
  colnames(design) <- unique_dpis

  fit <- lmFit(matrix_delta, design)
  fit_eb <- eBayes(fit, robust = TRUE)

  overall <- topTable(fit_eb, coef = seq_len(ncol(design)), number = Inf, sort.by = "none")
  overall$product <- rownames(overall)
  overall$threshold <- threshold
  overall <- overall[, c("threshold", "product", "AveExpr", "F", "P.Value", "adj.P.Val")]
  names(overall) <- c("threshold", "product", "ave_expr", "moderated_F", "p_value", "bh_fdr")

  if (write_detailed) {
    overall_file <- file.path(outdir, paste0("limma_overall_by_gene_", threshold, ".tsv"))
    write_tsv(overall, overall_file)
  }

  # Build dynamic contrasts based on unique DPIs
  if (length(unique_dpis) < 2) {
    message("Only 1 DPI group (", unique_dpis[1], ") found. Skipping contrast fitting for threshold ", threshold)
    # Return overall only
    return(list(overall = overall, contrasts = NULL, sample_info = sample_info))
  }

  contrast_exprs <- c()
  contrast_names <- c()

  # Successive contrasts (e.g. dpi2_vs_dpi1, dpi3_vs_dpi2...)
  for (i in 2:length(unique_dpis)) {
    c_name <- paste0(unique_dpis[i], "_vs_", unique_dpis[i-1])
    c_expr <- paste0(unique_dpis[i], " - ", unique_dpis[i-1])
    contrast_names <- c(contrast_names, c_name)
    contrast_exprs <- c(contrast_exprs, c_expr)
  }
  
  # Overall contrast (last vs first)
  if (length(unique_dpis) > 2) {
    c_name <- paste0(unique_dpis[length(unique_dpis)], "_vs_", unique_dpis[1])
    c_expr <- paste0(unique_dpis[length(unique_dpis)], " - ", unique_dpis[1])
    contrast_names <- c(contrast_names, c_name)
    contrast_exprs <- c(contrast_exprs, c_expr)
  }

  names(contrast_exprs) <- contrast_names
  
  # Parse contrasts in R environment
  contrasts_expr_parsed <- paste(contrast_names, "=", paste0("'", contrast_exprs, "'"), collapse = ", ")
  eval_str <- paste0("makeContrasts(", paste(paste0(contrast_names, " = ", contrast_exprs), collapse = ", "), ", levels = design)")
  contrasts_obj <- eval(parse(text = eval_str))

  fit2 <- contrasts.fit(fit, contrasts_obj)
  fit2 <- eBayes(fit2, robust = TRUE)

  contrast_tables <- list()
  for (coef_name in colnames(contrasts_obj)) {
    tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "none")
    tt$product <- rownames(tt)
    tt$threshold <- threshold
    tt$contrast <- coef_name
    tt <- tt[, c("threshold", "contrast", "product", "logFC", "t", "P.Value", "adj.P.Val", "B")]
    names(tt) <- c("threshold", "contrast", "product", "logFC", "moderated_t", "p_value", "bh_fdr", "B_stat")

    if (write_detailed) {
      contrast_file <- file.path(outdir, paste0("limma_contrast_", coef_name, "_", threshold, ".tsv"))
      write_tsv(tt, contrast_file)
    }
    contrast_tables[[coef_name]] <- tt
  }

  all_contrasts <- do.call(rbind, contrast_tables)
  if (write_detailed) {
    contrasts_file <- file.path(outdir, paste0("limma_contrasts_by_gene_", threshold, ".tsv"))
    write_tsv(all_contrasts, contrasts_file)
  }

  list(overall = overall, contrasts = all_contrasts, sample_info = sample_info)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  if (!file.exists(opts$delta_file)) {
    stop(
      paste0(
        "Missing delta file: ", opts$delta_file,
        "\nRun python3 SNPGenie/scripts/analyze_delta_selection.py first."
      )
    )
  }

  if (opts$manifest == "" || !file.exists(opts$manifest)) {
    stop("Missing or invalid samples_manifest.tsv path. Provide --manifest=<path>")
  }

  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

  # Read manifest
  manifest_df <- read.delim(opts$manifest, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

  df <- read.delim(opts$delta_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  df$delta_piN_minus_piS <- as.numeric(df$delta_piN_minus_piS)

  thresholds <- sort(unique(df$threshold))
  overall_all <- list()
  contrasts_all <- list()
  sample_info_all <- list()

  for (thr in thresholds) {
    df_thr <- df[df$threshold == thr, , drop = FALSE]
    out <- run_limma_for_threshold(df_thr, thr, opts$outdir, manifest_df, opts$write_detailed)
    
    if (is.null(out)) next
    
    overall_all[[thr]] <- out$overall
    if (!is.null(out$contrasts)) {
      contrasts_all[[thr]] <- out$contrasts
    }

    sample_info <- out$sample_info
    sample_info$threshold <- thr
    sample_info_all[[thr]] <- sample_info[, c("threshold", "sample", "dpi")]
  }

  if (length(overall_all) > 0) {
    overall_combined <- do.call(rbind, overall_all)
    write_tsv(overall_combined, file.path(opts$outdir, "limma_overall_by_gene_all_thresholds.tsv"))
    
    sample_info_combined <- do.call(rbind, sample_info_all)
    write_tsv(sample_info_combined, file.path(opts$outdir, "limma_sample_design.tsv"))
  }

  if (length(contrasts_all) > 0) {
    contrasts_combined <- do.call(rbind, contrasts_all)
    write_tsv(contrasts_combined, file.path(opts$outdir, "limma_contrasts_by_gene_all_thresholds.tsv"))
  } else {
    # Write empty contrasts file if none exists (e.g. single DPI dataset)
    empty_df <- data.frame(
      threshold = character(0),
      contrast = character(0),
      product = character(0),
      logFC = numeric(0),
      moderated_t = numeric(0),
      p_value = numeric(0),
      bh_fdr = numeric(0),
      B_stat = numeric(0)
    )
    write_tsv(empty_df, file.path(opts$outdir, "limma_contrasts_by_gene_all_thresholds.tsv"))
  }

  message("Wrote limma outputs to: ", opts$outdir)
}

main()
