#!/usr/bin/env Rscript

# Loads all packages in a way that allows exporting to child environments
packages <- c("ggplot2", "ggthemes")
for (p in packages) {
  library(p, character.only = TRUE)
}

######
# DATA PREP
######

setwd("C.albicans_in_vitro_pipeline/src/mhclust/")
source("clusterv2.R")

# Loads the trimmed mutant profiles with LFC>=1 in at least one condition across all but YPD
scores <- read.csv("../../data/output/summary_mod_t_test_LFC_normvar_LFCgeq1_noNA_noYPD.txt", sep = "\t", row.names = 1)
output_folder <- file.path("output_clusters_v2_LFCgeq1_noNA_noYPD")
if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }

# Sets parameters that determine the number of clusters for each run
n_layers <- 3
metric_thresholds <- c(0.01, 0.01, 0.01, 0.01)  # One more than the number of layers if we don't load existing clusters
metric <- "snr"
sort_metric <- "size"
sort_decreasing <- FALSE
member_threshold <- 3
#filter_file <- file.path(output_folder, "orig_4k_layer1_clusters_to_filter_full283_FDR5.txt")
filter_file <- NULL

# Sets whether or not to load existing clusters by giving a path to that file or generate new ones,
# by setting existing_clusters to NULL
existing_clusters <- NULL
#existing_clusters <- file.path(output_folder, "orig_layer1", "best_clusters.txt")

# Replaces NAs with row means
row_means <- rowMeans(scores, na.rm = TRUE)
for (i in 1:nrow(scores)) {
  na_ind <- which(is.na(scores[i,]))
  scores[i, na_ind] <- row_means[i]
}

######
# FIRST LAYER PREP
######

# Either loads in existing clusters or clusters qGI scores directly
first_clust <- NULL
layer_output_folder <- file.path(output_folder, "layer0")
if (!dir.exists(layer_output_folder)) { dir.create(layer_output_folder) }
if (is.null(existing_clusters)) {
  gi_cluster(scores, layer_output_folder,
             metric_threshold = metric_thresholds[1],
             metric = metric,
             sort_metric = sort_metric,
             sort_decreasing = sort_decreasing,
             pick_num = 1,
             member_threshold = member_threshold,
             verbose = TRUE)
  first_clust <- readLines(file.path(layer_output_folder, "best_clusters.txt"))
} else {
  first_clust <- readLines(existing_clusters)
}

# Manually filters the first layer of clusters if file specified
if (!is.null(filter_file)) {
  to_remove <- as.numeric(readLines(filter_file))
  to_keep <- 1:length(first_clust)
  to_keep <- to_keep[!(to_keep %in% to_remove)]
  first_clust <- first_clust[to_keep]
}

######
# MAIN SCRIPT
######

# Creates metagene dataframe for the first layer
layer_df <- data.frame(matrix(nrow = length(first_clust), ncol = ncol(scores)))
for (i in 1:length(first_clust)) {
  clust <- strsplit(first_clust[i], ";")[[1]]
  layer_df[i,] <- colMeans(scores[clust,])
}
rownames(layer_df) <- paste0("Clust0_", 1:nrow(layer_df))
colnames(layer_df) <- colnames(scores)
write.table(layer_df, file.path(layer_output_folder, "merged_scores.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Clusters each layer
previous_layer_clust <- first_clust
previous_layer_genes <- first_clust
previous_layer_df <- layer_df
for (layer in 1:n_layers) {

  # Makes output folder for the current layer
  layer_output_folder <- file.path(output_folder, paste0("layer", layer))
  if (!dir.exists(layer_output_folder)) { dir.create(layer_output_folder) }
  print(layer_df)
  # Clusters the current layer
  gi_cluster(layer_df, layer_output_folder,
             metric_threshold = metric_thresholds[layer+1],
             metric = metric,
             sort_metric = sort_metric,
             sort_decreasing = sort_decreasing,
             pick_num = 1,
             member_threshold = member_threshold,
             verbose = TRUE)

  # Loads clusters
  current_clust <- readLines(file.path(layer_output_folder, "best_clusters.txt"))

  # Creates metagene dataframe for the current layer
  layer_df <- data.frame(matrix(nrow = length(current_clust), ncol = ncol(scores)))
  for (i in 1:length(current_clust)) {
    clust <- strsplit(current_clust[i], ";")[[1]]
    layer_df[i,] <- colMeans(previous_layer_df[clust,])
  }
  clust_str <- paste0("Clust", layer, "_")
  rownames(layer_df) <- paste0(clust_str, 1:nrow(layer_df))
  colnames(layer_df) <- colnames(scores)
  write.table(layer_df, file.path(layer_output_folder, "merged_scores.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Maps clusters to gene names
  gene_list <- c()
  for (i in 1:length(current_clust)) {
    clust <- strsplit(current_clust[i], ";")[[1]]
    all_genes <- NULL
    for (id in clust) {
      clust_id <- as.numeric(strsplit(id, "_")[[1]][2])
      clust_genes <- strsplit(previous_layer_genes[clust_id], ";")[[1]]
      if(is.null(all_genes)) {
        all_genes <- clust_genes
      } else {
        all_genes <- c(all_genes, clust_genes)
      }
    }
    all_genes <- paste(sort(all_genes), collapse = ";")
    gene_list <- c(gene_list, all_genes)
  }
  writeLines(gene_list, file.path(layer_output_folder, "clust_genes.txt"))

  # Saves variables for next layer of clustering
  previous_layer_df <- layer_df
  previous_layer_clust <- current_clust
  previous_layer_genes <- gene_list
}

# Writes clustering data to zip
fname <- file.path(output_folder, paste0(basename(output_folder), ".zip"))
zip(fname, dir(output_folder, full.names = TRUE))
stop("Script finished")
