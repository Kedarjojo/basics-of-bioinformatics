source("config.R")

combined <- readRDS(
  file.path(OBJECTS_DIR, "combined_integrated.rds")
)

seurat_list <- SplitObject(combined, split.by = "dataset")

markers_list <- lapply(seurat_list, function(obj) {
  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- "seurat_clusters"
  FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
})

write.csv(markers_list[["GSE72056"]],
          file.path(TABLES_DIR, "markers_72056.csv"))

write.csv(markers_list[["GSE115978"]],
          file.path(TABLES_DIR, "markers_115978.csv"))