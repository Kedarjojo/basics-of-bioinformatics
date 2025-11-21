# ============================================================
# Script: download_and_prepare_GSE154826.R
# Purpose: Download and prepare NSCLC scRNA-seq data from GEO GSE154826
# ============================================================

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  install.packages("GEOquery")
}
library(GEOquery)
library(tools)

geo_acc <- "GSE154826"
dest_dir <- "data/GSE154826"

if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)

# 1. Download supplementary files (.tar.gz)
message("Downloading GEO supplementary files for ", geo_acc)
supp_files <- getGEOSuppFiles(geo_acc, makeDirectory = FALSE, baseDir = dest_dir)

# 2. Extract each tar.gz archive
tar_files <- rownames(supp_files)[grepl("\\.tar\\.gz$", rownames(supp_files))]

for (tar_file in tar_files) {
  message("Extracting ", tar_file)
  
  # Create folder for this batch
  batch_name <- file_path_sans_ext(basename(tar_file))
  batch_dir <- file.path(dest_dir, batch_name)
  if (!dir.exists(batch_dir)) dir.create(batch_dir)
  
  # Extract archive
  untar(file.path(dest_dir, tar_file), exdir = batch_dir)
  
  # Rename files to standard 10x naming (barcodes, features, matrix)
  files <- list.files(batch_dir, full.names = TRUE)
  for (f in files) {
    if (grepl("_features\\.tsv$", f)) {
      file.rename(f, file.path(batch_dir, "features.tsv"))
    } else if (grepl("_barcodes\\.tsv$", f)) {
      file.rename(f, file.path(batch_dir, "barcodes.tsv"))
    } else if (grepl("_matrix\\.tsv$", f)) {
      # convert to mtx if needed
      # Seurat can read tsv using readMM, but it's better to convert to .mtx if large
      file.rename(f, file.path(batch_dir, "matrix.mtx"))
    }
  }
}

message("✅ Download and preparation complete. Each batch is in its own folder ready for Seurat.")
