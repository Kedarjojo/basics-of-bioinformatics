source("scripts/config.R")

seurat_72056  <- readRDS(file.path(OBJECTS_DIR, "seurat_72056.rds"))
seurat_115978 <- readRDS(file.path(OBJECTS_DIR, "seurat_115978.rds"))

obj_list <- list(seurat_72056, seurat_115978)

obj_list <- lapply(obj_list, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  return(obj)
})

anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  dims = 1:30
)

combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:15)
combined <- FindNeighbors(combined, dims = 1:15)
combined <- FindClusters(combined)

saveRDS(combined,
        file.path(OBJECTS_DIR, "combined_integrated.rds"))