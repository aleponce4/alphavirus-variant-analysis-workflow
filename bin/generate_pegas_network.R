#!/usr/bin/env Rscript
# generate_pegas_network.R
# Generate publication-grade quasispecies haplotype Minimum Spanning Network using R pegas haploNet

suppressPackageStartupMessages({
  library(ape)
  library(pegas)
  library(igraph)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript generate_pegas_network.R <frequency_by_sample.csv> <output.png> [dataset_name]\n")
  quit(status = 1)
}

freq_csv <- args[1]
out_png  <- args[2]
dataset  <- ifelse(length(args) >= 3, args[3], "Haplotype Network")

if (!file.exists(freq_csv)) {
  cat("Notice: Missing frequency input CSV for pegas network plot. Skipping.\n")
  quit(status = 0)
}

# 1. Read frequency table
freq_df <- read.csv(freq_csv, stringsAsFactors = FALSE)

# Filter out massive deletion artifacts (defective interfering particles) that break the scale
freq_df <- freq_df[is.na(freq_df$n_mutations) | freq_df$n_mutations <= 15, ]

# 2. Filter to qualifying haplotypes (frequency >= 0.01 in at least one sample)
qual_haps <- unique(freq_df$haplotype_id[freq_df$frequency >= 0.01])
qual_haps <- sort(qual_haps)

if (length(qual_haps) < 2) {
  cat("Notice: Fewer than 2 qualifying haplotypes found. Skipping pegas network plot.\n")
  quit(status = 0)
}

sub_df <- freq_df[freq_df$haplotype_id %in% qual_haps, ]

# 3. Extract mutation signatures per haplotype and handle Reference
hap_muts <- list()
ref_hap <- NULL

# Identify if any existing haplotype IS the reference (0 mutations)
for (h in qual_haps) {
  mut_str <- sub_df$mutations[sub_df$haplotype_id == h][1]
  if (is.na(mut_str) || mut_str == "REF" || mut_str == "") {
    ref_hap <- h
    break
  }
}

# Rename or Inject Reference
if (!is.null(ref_hap)) {
  qual_haps[qual_haps == ref_hap] <- "Reference"
  sub_df$haplotype_id[sub_df$haplotype_id == ref_hap] <- "Reference"
} else {
  qual_haps <- c("Reference", qual_haps)
}

# Build mutation list
for (h in qual_haps) {
  if (h == "Reference") {
    hap_muts[[h]] <- character(0)
  } else {
    mut_str <- sub_df$mutations[sub_df$haplotype_id == h][1]
    hap_muts[[h]] <- unlist(strsplit(mut_str, ";|,"))
  }
}

# 4. Compute SNV mutation count distance matrix
n_nodes <- length(qual_haps)
d_mat <- matrix(0, nrow = n_nodes, ncol = n_nodes, dimnames = list(qual_haps, qual_haps))

for (i in 1:n_nodes) {
  for (j in 1:n_nodes) {
    if (i != j) {
      m1 <- hap_muts[[qual_haps[i]]]
      m2 <- hap_muts[[qual_haps[j]]]
      dist_val <- length(setdiff(m1, m2)) + length(setdiff(m2, m1))
      d_mat[i, j] <- dist_val
    }
  }
}

d_dist <- as.dist(d_mat)

# 5. Build pegas haploNet network via reticulate minimum spanning tree (PopART style)
net <- pegas::rmst(d_dist, quiet = TRUE)

# 6. Skip igraph force-directed layout to allow pegas to use its native 
# PopART-style distance-based geometric layout.
# ig <- pegas::as.igraph.haploNet(net)
# xy_mat <- igraph::layout_with_fr(ig)
# xy_coords <- data.frame(x = xy_mat[, 1], y = xy_mat[, 2])

# 7. Okabe-Ito color palette
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

group_col <- "treatment"
group_title <- "Group"

if ("dpi" %in% colnames(sub_df) && length(unique(sub_df$dpi[!is.na(sub_df$dpi) & sub_df$dpi != ""])) > 1) {
  group_col <- "dpi"
  group_title <- "DPI"
} else if ("day" %in% colnames(sub_df) && length(unique(sub_df$day[!is.na(sub_df$day) & sub_df$day != ""])) > 1) {
  group_col <- "day"
  group_title <- "Day"
}

groups <- sort(unique(sub_df[[group_col]]))
groups <- groups[groups != "" & !is.na(groups)]
if (length(groups) == 0) groups <- "group"
group_colors <- setNames(okabe_ito[1:length(groups)], groups)

# 8. Build pie matrix for group composition per haplotype
net_labels <- labels(net)
pie_matrix <- matrix(0, nrow = length(net_labels), ncol = length(groups),
                     dimnames = list(net_labels, groups))

for (i in seq_along(net_labels)) {
  node_name <- net_labels[i]
  if (node_name == "Reference" && is.null(ref_hap)) {
    # Injected reference (not in samples) gets tiny artificial size
    for (g in groups) pie_matrix[i, g] <- 0
  } else {
    node_df <- sub_df[sub_df$haplotype_id == node_name, ]
    for (g in groups) {
      g_freqs <- node_df$frequency[node_df[[group_col]] == g]
      pie_matrix[i, g] <- if (length(g_freqs) > 0) sum(g_freqs, na.rm = TRUE) else 0
    }
  }
}

pie_matrix[is.na(pie_matrix)] <- 0
row_sums <- rowSums(pie_matrix)

# If any node has 0 abundance (e.g. injected reference), give it a tiny dummy value
zero_rows <- which(row_sums == 0)
if (length(zero_rows) > 0) {
  pie_matrix[zero_rows, 1] <- 0.02
  row_sums[zero_rows] <- 0.02
}

# Store the raw sums for size scaling
sizes <- row_sums
sizes[sizes == 0] <- 0.02 # ensure injected reference has minimal size
# Use strictly proportional area scaling (no flat intercept)
sizes_scaled <- sqrt(sizes) * 3.5

# Normalize pie_matrix so rows sum to 1 to prevent "incomplete circles" in pegas
for (i in 1:nrow(pie_matrix)) {
  if (row_sums[i] > 0) {
    pie_matrix[i, ] <- pie_matrix[i, ] / row_sums[i]
  }
}

# Render 5-inch wide publication figure using Arial font
png(out_png, width = 5, height = 4.5, units = "in", res = 300, family = "Arial")
par(mar = c(1.5, 1.5, 2.5, 1.5), family = "Arial", cex.main = 0.9)

plot(
  net,
  size = sizes_scaled,
  pie = pie_matrix,
  bg = group_colors,
  # xy = xy_coords, # Removed to use PopART geometric layout
  fast = FALSE,
  labels = TRUE,
  font = 2,
  cex = 0.5,
  threshold = 0,
  show.mutation = 1, # PopART hatch marks
  scale.ratio = 1.0, # PopART proportional segments
  main = paste0(dataset, " — Quasispecies Temporal Haplotype Network (PopART style)")
)

# 1. Group / DPI Legend
leg1 <- legend("topleft", legend = groups, fill = group_colors, bty = "n", cex = 0.7, title = group_title)

# Determine dynamic legend values based on max abundance
max_abund <- max(sizes)
if (max_abund <= 1.0) {
  leg_vals <- c(0.1, 0.5, 1.0)
} else if (max_abund <= 3.0) {
  leg_vals <- c(0.1, 1.0, ceiling(max_abund))
} else {
  leg_vals <- c(1.0, floor(max_abund / 2), ceiling(max_abund))
}
leg_vals <- unique(leg_vals)

# 2. Total Abundance Node Size Legend
# Placed dynamically below the first legend, with increased vertical spacing (y.intersp)
legend(
  x = leg1$rect$left,
  y = leg1$rect$top - leg1$rect$h - (diff(par("usr")[3:4]) * 0.04),
  legend = as.character(leg_vals),
  pt.cex = sqrt(leg_vals) * 3.5,
  pch = 21,
  col = "black",
  pt.bg = "white",
  bty = "n",
  title = "Total Abundance",
  cex = 0.7,
  y.intersp = 2.5,
  x.intersp = 2.5
)

dev.off()
cat("Wrote pegas haplotype network plot:", out_png, "\n")
