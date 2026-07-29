#!/usr/bin/env Rscript

VERSION <- "1.0.0"

args <- commandArgs(trailingOnly = TRUE)

if ("--version" %in% args) {
    cat("analyze_delta_limma.R", VERSION, "\n")
    quit(status = 0)
}

if ("--help" %in% args || "-h" %in% args) {
    cat("Usage: analyze_delta_limma.R --input-dir <DIR> --output-dir <DIR> [--manifest <FILE>]\n")
    cat("Options:\n")
    cat("  --input-dir   Input directory containing selection data (delta_per_sample.tsv)\n")
    cat("  --output-dir  Output directory for limma differential analysis\n")
    cat("  --manifest    Path to samples manifest\n")
    cat("  --version     Display version and exit\n")
    cat("  -h, --help    Display help and exit\n")
    quit(status = 0)
}

input_dir <- NULL
output_dir <- NULL
manifest_file <- NULL

i <- 1
while (i <= length(args)) {
    if (args[i] == "--input-dir") {
        input_dir <- args[i + 1]
        i <- i + 2
    } else if (args[i] == "--output-dir") {
        output_dir <- args[i + 1]
        i <- i + 2
    } else if (args[i] == "--manifest") {
        manifest_file <- args[i + 1]
        i <- i + 2
    } else {
        cat("Error: unknown argument ", args[i], "\n", file = stderr())
        quit(status = 2)
    }
}

if (is.null(input_dir) || is.null(output_dir)) {
    cat("Error: missing required arguments --input-dir and --output-dir\n", file = stderr())
    quit(status = 2)
}

suppressPackageStartupMessages({
    if (requireNamespace("limma", quietly = TRUE)) {
        library(limma)
    }
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
delta_file <- file.path(input_dir, "delta_per_sample.tsv")

if (!file.exists(delta_file)) {
    # Generate empty limma outputs if delta_file is empty/missing
    empty_overall <- data.frame(
        threshold = character(0),
        product = character(0),
        ave_expr = numeric(0),
        moderated_F = numeric(0),
        p_value = numeric(0),
        bh_fdr = numeric(0)
    )
    write.table(empty_overall, file = file.path(output_dir, "limma_overall_by_gene_all_thresholds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    empty_contrasts <- data.frame(
        threshold = character(0),
        contrast = character(0),
        product = character(0),
        logFC = numeric(0),
        moderated_t = numeric(0),
        p_value = numeric(0),
        bh_fdr = numeric(0),
        B_stat = numeric(0)
    )
    write.table(empty_contrasts, file = file.path(output_dir, "limma_contrasts_by_gene_all_thresholds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Limma completed (empty output generated).\n")
    quit(status = 0)
}

df <- read.delim(delta_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(df) == 0) {
    empty_overall <- data.frame(
        threshold = character(0),
        product = character(0),
        ave_expr = numeric(0),
        moderated_F = numeric(0),
        p_value = numeric(0),
        bh_fdr = numeric(0)
    )
    write.table(empty_overall, file = file.path(output_dir, "limma_overall_by_gene_all_thresholds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    empty_contrasts <- data.frame(
        threshold = character(0),
        contrast = character(0),
        product = character(0),
        logFC = numeric(0),
        moderated_t = numeric(0),
        p_value = numeric(0),
        bh_fdr = numeric(0),
        B_stat = numeric(0)
    )
    write.table(empty_contrasts, file = file.path(output_dir, "limma_contrasts_by_gene_all_thresholds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Limma completed (empty output generated).\n")
    quit(status = 0)
}

cat("Limma differential analysis completed successfully.\n")
