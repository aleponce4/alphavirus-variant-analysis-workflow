#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
})

parse_args <- function(args) {
  opts <- list(
    delta_file = "SNPGenie/analysis/delta_selection/delta_per_sample.tsv",
    outdir = "SNPGenie/analysis/delta_selection",
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
    } else if (arg == "--write-detailed") {
      opts$write_detailed <- TRUE
    }
  }
  opts
}

extract_dpi <- function(sample_name) {
  m <- regexec("INH_([0-9]+)_DPI_", sample_name)
  reg <- regmatches(sample_name, m)
  if (length(reg[[1]]) >= 2) {
    return(paste0("dpi", reg[[1]][2]))
  }
  return(NA_character_)
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

run_limma_for_threshold <- function(df_thr, threshold, outdir, write_detailed = FALSE) {
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

  sample_info <- data.frame(
    sample = colnames(matrix_delta),
    dpi = vapply(colnames(matrix_delta), extract_dpi, character(1)),
    stringsAsFactors = FALSE
  )

  sample_info <- sample_info[sample_info$dpi %in% c("dpi1", "dpi3", "dpi5"), , drop = FALSE]
  sample_info$dpi <- factor(sample_info$dpi, levels = c("dpi1", "dpi3", "dpi5"))

  matrix_delta <- matrix_delta[, sample_info$sample, drop = FALSE]

  design <- model.matrix(~ 0 + dpi, data = sample_info)
  colnames(design) <- levels(sample_info$dpi)

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

  contrasts <- makeContrasts(
    dpi3_vs_dpi1 = dpi3 - dpi1,
    dpi5_vs_dpi3 = dpi5 - dpi3,
    dpi5_vs_dpi1 = dpi5 - dpi1,
    levels = design
  )

  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2, robust = TRUE)

  contrast_tables <- list()
  for (coef_name in colnames(contrasts)) {
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

  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

  df <- read.delim(opts$delta_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  df <- df[!grepl("pilot", df$sample, ignore.case = TRUE), , drop = FALSE]
  df$delta_piN_minus_piS <- as.numeric(df$delta_piN_minus_piS)

  thresholds <- sort(unique(df$threshold))
  overall_all <- list()
  contrasts_all <- list()
  sample_info_all <- list()

  for (thr in thresholds) {
    df_thr <- df[df$threshold == thr, , drop = FALSE]
    out <- run_limma_for_threshold(df_thr, thr, opts$outdir, opts$write_detailed)
    overall_all[[thr]] <- out$overall
    contrasts_all[[thr]] <- out$contrasts

    sample_info <- out$sample_info
    sample_info$threshold <- thr
    sample_info_all[[thr]] <- sample_info[, c("threshold", "sample", "dpi")]
  }

  overall_combined <- do.call(rbind, overall_all)
  contrasts_combined <- do.call(rbind, contrasts_all)
  sample_info_combined <- do.call(rbind, sample_info_all)

  write_tsv(overall_combined, file.path(opts$outdir, "limma_overall_by_gene_all_thresholds.tsv"))
  write_tsv(contrasts_combined, file.path(opts$outdir, "limma_contrasts_by_gene_all_thresholds.tsv"))
  write_tsv(sample_info_combined, file.path(opts$outdir, "limma_sample_design.tsv"))

  message("Wrote limma outputs to: ", opts$outdir)
}

main()
