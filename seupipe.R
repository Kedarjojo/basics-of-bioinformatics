library(Seurat) 
library(Matrix) 
library(dplyr) 
library(ggplot2)

load_scRNA_data <- function(directory) {
  barcode_file <- file.path(directory, "barcodes.tsv")
  feature_file <- file.path(directory, "features.tsv")
  matrix_file <- file.path(directory, "matrix.mtx")
  
  # Load data
  sparse_matrix <- readMM(matrix_file)
  barcodes <- readLines(barcode_file)
  features <- read.table(feature_file, header = FALSE, stringsAsFactors = FALSE)
  
  # Check for duplicate features
  if (anyDuplicated(features$V2)) {
    warning("Duplicate feature names found. Renaming to ensure uniqueness.")
    features$V2 <- make.unique(features$V2)  # Ensure unique feature names
  }
  
  # Set row and column names
  rownames(sparse_matrix) <- features$V2
  colnames(sparse_matrix) <- barcodes
  
  # Create Seurat object
  seurat_object <- CreateSeuratObject(counts = sparse_matrix, project = basename(directory))
  return(seurat_object)
}


directories <- "C:/Users/kjoshi/Downloads/filtered_feature_bc_matrix"
seurat_list <- lapply(directories, load_scRNA_data)


combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1]) 
combined_seurat <- NormalizeData(combined_seurat) 
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000) 
combined_seurat <- ScaleData(combined_seurat) 
combined_seurat <- RunPCA(combined_seurat) 
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:10) 
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5) 
combined_seurat <- RunUMAP(combined_seurat, dims = 1:10) 


DefaultAssay(combined_seurat)
features <- Features(combined_seurat@assays$RNA)



DimPlot(combined_seurat, reduction = "umap", group.by = "seurat_clusters") + ggtitle("UMAP Clustering")
