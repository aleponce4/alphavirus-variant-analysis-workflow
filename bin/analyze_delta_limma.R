#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    if (!requireNamespace("optparse", quietly = TRUE)) {
        # Fallback for stub execution if optparse not yet installed in host environment
    }
})

VERSION <- "1.0.0"

args <- commandArgs(trailingOnly = TRUE)

if ("--version" %in% args) {
    cat("analyze_delta_limma.R", VERSION, "\n")
    quit(status = 0)
}

if ("--help" %in% args || "-h" %in% args) {
    cat("Usage: analyze_delta_limma.R --input-dir <DIR> --output-dir <DIR>\n")
    cat("Options:\n")
    cat("  --input-dir   Input directory containing selection data\n")
    cat("  --output-dir  Output directory for limma differential analysis\n")
    cat("  --version     Display version and exit\n")
    cat("  -h, --help    Display help and exit\n")
    quit(status = 0)
}

# Check for unexpected arguments or missing required arguments
input_dir <- NULL
output_dir <- NULL

i <- 1
while (i <= length(args)) {
    if (args[i] == "--input-dir") {
        input_dir <- args[i + 1]
        i <- i + 2
    } else if (args[i] == "--output-dir") {
        output_dir <- args[i + 1]
        i <- i + 2
    } else {
        cat("Error: unknown argument ", args[i], "\n", file = stderr())
        quit(status = 2)
    }
}

cat("analyze_delta_limma.R stub executed successfully.\n")
