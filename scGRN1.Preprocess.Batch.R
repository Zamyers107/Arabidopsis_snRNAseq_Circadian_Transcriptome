# ------------------------------
# Libraries and Setup
# ------------------------------
.libPaths("/home/greenhamlab/bin/R-4.4.0/")
library(readxl)
library(dplyr)
library(tidyverse)
library(Seurat)
library(SCopeLoomR)
library(tools)

# ------------------------------
# Configuration
# ------------------------------
seuratObj_path <- "/home/greenhamlab/snTC/20250127.snTC.GRNs.Subclustered/20241217.snTC.DROP24B.Subclustered.c8c12.RDS"
cycling_folder <- "/home/greenhamlab/snTC/20250227.snTC.ClusterSpecificGRNs/0.05p_sorted_JTK_Data/"
group_col <- "sub.cluster"
today <- format(Sys.Date(), "%Y%m%d")

# Clusters to process
all_clusters <- c("0", "1", "5", "2", "11", "3", "4", "12_2", "10", "9", "14", "6", 
                  "12_0", "7", "8_0", "12_1", "13", "8_1", "15") # Excludes "16"

# Base directory for output
base_outdir <- "/home/greenhamlab/snTC/20240411.GRNs.ParallelShuffled"

# ------------------------------
# Helper Functions
# ------------------------------
as_matrix <- function(mat){
  tmp <- matrix(data = 0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i + 1
  col_pos <- findInterval(seq(mat@x)-1, mat@p[-1]) + 1
  val <- mat@x
  for (i in seq_along(val)) {
    tmp[row_pos[i], col_pos[i]] <- val[i]
  }
  rownames(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

add_cellAnnotation <- function(loom, cellAnnotation){
  cellAnnotation <- data.frame(cellAnnotation)
  cellAnnotation <- cellAnnotation[, !colnames(cellAnnotation) %in% c("nGene", "nUMI"), drop = FALSE]
  cellAnnotation <- cellAnnotation[get_cell_ids(loom), , drop = FALSE]
  for (cn in colnames(cellAnnotation)) {
    add_col_attr(loom = loom, key = cn, value = cellAnnotation[, cn])
  }
  invisible(loom)
}

# ------------------------------
# Load Data
# ------------------------------
SeuratObj <- readRDS(seuratObj_path)
file_names <- list.files(cycling_folder, pattern = "*.xlsx", full.names = TRUE)
cluster_data <- lapply(file_names, read_excel)
names(cluster_data) <- gsub("_.*$", "", basename(file_names))

# ------------------------------
# Main Loop
# ------------------------------
for (cluster_id in all_clusters) {
  message("Processing cluster: ", cluster_id)
  outdir <- file.path(base_outdir, paste0("GRNs.c", cluster_id, ".ClusterOnlyCycling.c", cluster_id, "OnlyCells"))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset cluster from the original Seurat object
  SeuratCluster <- subset(SeuratObj, sub.cluster == cluster_id)
  cellmeta <- SeuratCluster@meta.data
  
  # Get the expression matrix from the cluster (sparse matrix)
  mat <- GetAssayData(SeuratCluster, layer = "data")
  
  # ------------------------------
  # Filter for cluster-specific cycling genes
  # ------------------------------
  cluster_key <- paste0("c", cluster_id)
  if (!cluster_key %in% names(cluster_data)) {
    warning("No cycling gene data found for cluster: ", cluster_id)
    next
  }
  unique_cycid <- distinct(select(cluster_data[[cluster_key]], CycID))
  colnames(unique_cycid) <- c("gid")
  
  # Determine the genes to keep
  genes_to_keep <- intersect(unique_cycid$gid, rownames(mat))
  
  # ------------------------------
  # Create and Save Filtered Seurat Object
  # ------------------------------
  # Subset the Seurat object so that it is identical to the expression matrix used for the loom file.
  SeuratCluster_filtered <- subset(SeuratCluster, features = genes_to_keep)
  saveRDS(SeuratCluster_filtered, file = file.path(outdir, "Seurat_filtered.rds"))
  
  # ------------------------------
  # Generate the Expression Matrix and Save Looms
  # ------------------------------
  # Filter matrix to include only the selected genes
  mat_filtered <- mat[genes_to_keep, , drop = FALSE]
  
  # Remove the original SeuratCluster object (optional, for memory)
  rm(SeuratCluster)
  
  # Save cell metadata for additional processing
  saveRDS(cellmeta, file.path(outdir, "cellmeta.rds"))
  
  # Convert the filtered matrix (typically in a sparse format) to a dense matrix
  mat_filtered <- as_matrix(mat_filtered)
  
  # Save the original loom file
  loom <- build_loom(file.path(outdir, "exprMat.loom"), dgem = mat_filtered)
  loom <- add_cellAnnotation(loom, cellmeta)
  close_loom(loom)
  
  # ------------------------------
  # Create 10 Shuffled Looms
  # ------------------------------
  for (i in 1:10) {
    message("  Generating shuffled loom ", i, " for cluster ", cluster_id)
    # Shuffle each cell's values independently (each column)
    shuffled_mat <- apply(mat_filtered, 2, function(x) sample(x))
    # Reassign the original rownames so that gene identifiers are preserved
    rownames(shuffled_mat) <- rownames(mat_filtered)
    
    # Save shuffled loom file with a unique filename for each shuffle iteration
    shuffle_filename <- file.path(outdir, paste0("exprMat.shuffle", i, ".loom"))
    loom_shuffled <- build_loom(shuffle_filename, dgem = shuffled_mat)
    loom_shuffled <- add_cellAnnotation(loom_shuffled, cellmeta)
    close_loom(loom_shuffled)
  }
}

##### Capturing missed clusters #####

# ------------------------------
# Configuration
# ------------------------------
seuratObj_path <- "/home/greenhamlab/snTC/20250127.snTC.GRNs.Subclustered/20241217.snTC.DROP24B.Subclustered.c8c12.RDS"
cycling_folder <- "/home/greenhamlab/snTC/20250227.snTC.ClusterSpecificGRNs/0.05p_sorted_JTK_Data/"
group_col      <- "sub.cluster"
base_outdir    <- "/home/greenhamlab/snTC/20240411.GRNs.ParallelShuffled"

# Only re‑process these sub‑clusters:
missed_clusters <- c("8_0", "8_1", "12_0", "12_1", "12_2")

# ------------------------------
# Helper Functions
# ------------------------------
as_matrix <- function(mat) {
  tmp <- matrix(0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i + 1
  col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1]) + 1
  for (i in seq_along(mat@x)) {
    tmp[row_pos[i], col_pos[i]] <- mat@x[i]
  }
  rownames(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  tmp
}

add_cellAnnotation <- function(loom, cellAnnotation) {
  cellAnnotation <- as.data.frame(cellAnnotation)[, setdiff(colnames(cellAnnotation), c("nGene","nUMI")), drop = FALSE]
  cellAnnotation <- cellAnnotation[get_cell_ids(loom), , drop = FALSE]
  for (cn in colnames(cellAnnotation)) {
    add_col_attr(loom = loom, key = cn, value = cellAnnotation[[cn]])
  }
  invisible(loom)
}

# ------------------------------
# Load Seurat Object
# ------------------------------
message("Loading Seurat object from: ", seuratObj_path)
SeuratObj <- readRDS(seuratObj_path)

# ------------------------------
# Load & Name Cycling Data (fixed)
# ------------------------------
file_names   <- list.files(cycling_folder, pattern = "\\.xlsx$", full.names = TRUE)
cluster_data <- lapply(file_names, read_excel)

# strip off “_0.05” suffix (or any trailing “_<something>”) so names become "c8_0", "c12_2", etc.
raw_names   <- file_path_sans_ext(basename(file_names))
cluster_ids <- sub("^([cC]\\d+(?:_\\d+)?)_.*$", "\\1", raw_names)

names(cluster_data) <- cluster_ids
message("Found cycling data for clusters: ", paste(unique(cluster_ids), collapse = ", "))

# ------------------------------
# Process Only Missed Sub‑Clusters
# ------------------------------
for (cluster_id in missed_clusters) {
  message(">> Processing sub‑cluster: ", cluster_id)
  
  # build key ("c8_0", etc.)
  cluster_key <- paste0("c", cluster_id)
  
  if (! cluster_key %in% names(cluster_data)) {
    warning("   No cycling data for ", cluster_key, "; skipping.")
    next
  }
  
  # subset Seurat by sub.cluster
  SeuratCluster <- subset(SeuratObj, subset = !!sym(group_col) == cluster_id)
  cellmeta      <- SeuratCluster@meta.data
  mat           <- GetAssayData(SeuratCluster, layer = "data")
  
  # get cycling gene IDs
  cycids         <- distinct(select(cluster_data[[cluster_key]], CycID))
  colnames(cycids) <- "gid"
  genes_to_keep  <- intersect(cycids$gid, rownames(mat))
  
  # prepare output directory
  outdir <- file.path(base_outdir,
                      paste0("GRNs.c", cluster_id, ".ClusterOnlyCycling.c", cluster_id, "OnlyCells"))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # save filtered Seurat object
  Seurat_filt <- subset(SeuratCluster, features = genes_to_keep)
  saveRDS(Seurat_filt, file = file.path(outdir, "Seurat_filtered.rds"))
  
  # convert and save expression matrix + metadata
  mat_filt <- as_matrix(mat[genes_to_keep, , drop = FALSE])
  saveRDS(cellmeta, file = file.path(outdir, "cellmeta.rds"))
  
  # build & save the original loom
  loom_orig <- build_loom(file.path(outdir, "exprMat.loom"), dgem = mat_filt)
  loom_orig <- add_cellAnnotation(loom_orig, cellmeta)
  close_loom(loom_orig)
  
  # create 10 shuffled looms
  for (i in seq_len(10)) {
    message("   Shuffling iteration: ", i)
    shuffled <- apply(mat_filt, 2, sample)
    rownames(shuffled) <- rownames(mat_filt)
    
    loom_s <- build_loom(file.path(outdir, paste0("exprMat.shuffle", i, ".loom")),
                         dgem = shuffled)
    loom_s <- add_cellAnnotation(loom_s, cellmeta)
    close_loom(loom_s)
  }
}

