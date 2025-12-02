library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(SingleR)
library(celldex)   
library(monocle3)
library(patchwork)


# ============================================================
# NSCLC scRNA-seq Analysis Pipeline
# Using GEO dataset GSE154826
# ============================================================


# -----------------------------
# Helper function: Load individual 10x-style datasets
# -----------------------------
load_scRNA_data <- function(directory) {
  # Find files in the folder
  barcode_file <- list.files(directory, pattern = "barcodes.tsv$", full.names = TRUE)
  feature_file <- list.files(directory, pattern = "features.tsv$", full.names = TRUE)
  matrix_file <- list.files(directory, pattern = "matrix.mtx$", full.names = TRUE)
  
  if(length(barcode_file) == 0 | length(feature_file) == 0 | length(matrix_file) == 0) {
    stop(paste("Missing files in", directory))
  }
  
  # Read data
  sparse_matrix <- readMM(matrix_file)
  barcodes <- readLines(barcode_file)
  features <- read.table(feature_file, header = FALSE, stringsAsFactors = FALSE)
  
  # Ensure unique feature names
  if (anyDuplicated(features$V2)) {
    features$V2 <- make.unique(features$V2)
  }
  
  rownames(sparse_matrix) <- features$V2
  colnames(sparse_matrix) <- barcodes
  
  seurat_obj <- CreateSeuratObject(counts = sparse_matrix, project = basename(directory))
  return(seurat_obj)
}


# -----------------------------
# Step 1: Load all datasets
# -----------------------------
# Parent directory containing all batch folders
data_dir <- "/Users/kedarjoshi/Downloads/Traj_Analysis"

# Get all batch folders (non-recursive)
folders <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)

# Filter to only include the GSE154826 batch folders
folders <- folders[grepl("GSE154826_amp_batch_ID_", folders)]

# Load each folder into a Seurat object
seurat_list <- lapply(folders, load_scRNA_data)
saveRDS(seurat_list, file = "seurat_list_NSCLC_GSE154826.rds")


# -----------------------------
# Step 2: QC and Filtering
# -----------------------------
message("Performing QC and filtering")
for (i in 1:length(seurat_list)) {
  obj <- seurat_list[[i]]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
  
  # Filter low quality cells
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &
                  percent.mt < 15 & percent.ribo < 50)
  seurat_list[[i]] <- obj
  message("Dataset ", i, " filtered: ", ncol(obj), " cells retained")
}


# -----------------------------
# Step 3: Normalization, HVG, Scaling, PCA
# -----------------------------
message("Normalizing and identifying variable features")
for (i in 1:length(seurat_list)) {
  obj <- seurat_list[[i]]
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  seurat_list[[i]] <- obj
}


# -----------------------------
# Step 4: Doublet Detection
# -----------------------------
message("Running DoubletFinder")
for (i in 1:length(seurat_list)) {
  obj <- seurat_list[[i]]
  sweep.res.list <- paramSweep_v3(obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  nExp <- round(0.075 * ncol(obj)) # assume 7.5% doublets
  obj <- doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp)
  seurat_list[[i]] <- subset(obj, subset = DF.classifications_0.25_0.01_XXX == "Singlet") # adjust column name if needed
}


# -----------------------------
# Step 5: Integration
# -----------------------------
message("Integrating datasets")
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)


# -----------------------------
# Step 6: Scaling, PCA, Clustering, UMAP
# -----------------------------
message("Scaling and clustering integrated data")
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:20)


# -----------------------------
# Step 7: Cell Type Annotation (SingleR)
# -----------------------------
message("Annotating cell types using SingleR")
human_ref <- HumanPrimaryCellAtlasData()
singleR_res <- SingleR(test = combined@assays$RNA@data, ref = human_ref, labels = human_ref$label.main)
combined$SingleR_label <- singleR_res$labels


# -----------------------------
# Step 8: Marker Gene Detection
# -----------------------------
message("Finding cluster marker genes")
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# -----------------------------
# Step 9: Trajectory Analysis (Monocle3)
# -----------------------------
message("Running trajectory analysis")
cds <- as.cell_data_set(combined)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
# Optional: define root cells
# root_cells <- colnames(combined)[combined$SingleR_label=="Epithelial"]
# cds <- order_cells(cds, root_cells=root_cells)


# -----------------------------
# Step 10: Visualizations
# -----------------------------
message("Generating plots")
p1 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters") + ggtitle("UMAP Clustering")
p2 <- DimPlot(combined, reduction = "umap", group.by = "SingleR_label") + ggtitle("Cell Type Annotation")
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = TRUE, label_leaves = TRUE, label_branch_points = TRUE)

# Display plots
p1 + p2

