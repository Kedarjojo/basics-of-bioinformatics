library(Seurat)
library(dplyr)

# ----------------------------
# Step 1: Organize dataset paths
# ----------------------------
dataset_dirs <- list.files("path_to_datasets/", full.names = TRUE)
dataset_names <- basename(dataset_dirs)

# ----------------------------
# Step 2: Load and preprocess datasets individually
# ----------------------------
seurat_list <- list()
for (i in seq_along(dataset_dirs)) {
  data <- Read10X(data.dir = dataset_dirs[i])
  seurat_obj <- CreateSeuratObject(counts = data, project = dataset_names[i])
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  seurat_list[[dataset_names[i]]] <- seurat_obj
}

# ----------------------------
# Step 3: Batch integration
# ----------------------------
# Define batch size
batch_size <- 2
batches <- split(seurat_list, ceiling(seq_along(seurat_list)/batch_size))

integrated_batches <- list()

for (i in seq_along(batches)) {
  # Find anchors for the batch
  anchors <- FindIntegrationAnchors(object.list = batches[[i]], normalization.method = "SCT")
  # Integrate the batch
  integrated_batches[[i]] <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
}

# ----------------------------
# Step 4: Merge all integrated batches
# ----------------------------
# The merge function combines Seurat objects
seurat_combined <- Reduce(function(x, y) merge(x, y), integrated_batches)

# ----------------------------
# Step 5: Downstream analysis
# ----------------------------
seurat_combined <- RunPCA(seurat_combined, verbose = FALSE)
seurat_combined <- RunUMAP(seurat_combined, dims = 1:10)
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:10)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)
