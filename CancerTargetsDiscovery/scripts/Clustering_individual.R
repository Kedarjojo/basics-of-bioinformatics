
# Load config and Seurat objects
source("scripts/config.R")

# Assign markers per dataset
markers_72056  <- read.csv(
  file.path(TABLES_DIR, "markers_72056.csv"))
markers_115978 <- read.csv(
  file.path(TABLES_DIR, "markers_115978.csv"))

# Filter top markers
filtered_markers_72056 <- markers_72056 %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.25)

top_markers_72056 <- filtered_markers_72056 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

filtered_markers_115978 <- markers_115978 %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.25)

top_markers_115978 <- filtered_markers_115978 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

# Save filtered markers
write.csv(filtered_markers_72056,
          file.path(TABLES_DIR, "filtered_markers_72056.csv"),
          row.names = FALSE)

write.csv(filtered_markers_115978,
          file.path(TABLES_DIR, "filtered_markers_115978.csv"),
          row.names = FALSE)

# Manual cluster-to-celltype mapping
cluster2celltype_72056 <- c(
  "0"="T_cells_helper", "1"="T_cells_helper", "2"="B_cells",
  "3"="T_cells_cytotoxic", "4"="Tumor", "5"="T_cells_regulatory",
  "6"="Tumor_proliferating", "7"="Monocytes_Macrophages", "8"="Fibroblasts",
  "9"="Endothelial", "10"="NK_cells", "11"="Proliferating",
  "12"="Tumor", "13"="Fibroblasts", "14"="Tumor", "15"="B_cells_plasma",
  "16"="B_cells_naive", "17"="Fibroblasts", "18"="Tumor",
  "19"="Endothelial", "20"="Endothelial", "22"="Oligodendrocytes",
  "23"="Tumor", "24"="Fibroblasts"
)

cluster2celltype_115978 <- c(
  "0"="T_cells_helper", "1"="T_cells_helper", "2"="B_cells",
  "3"="T_cells_cytotoxic", "4"="Tumor", "5"="T_cells_regulatory",
  "6"="Tumor_proliferating", "7"="Monocytes_Macrophages", "8"="Fibroblasts",
  "9"="Endothelial", "10"="NK_cells", "11"="Proliferating",
  "12"="Tumor", "13"="Fibroblasts", "14"="Tumor", "15"="B_cells_plasma",
  "16"="B_cells_naive", "17"="Fibroblasts", "18"="Tumor",
  "19"="Endothelial", "20"="Endothelial", "21"="Tumor",
  "22"="Fibroblasts", "23"="Tumor", "24"="Fibroblasts"
)

map_clusters_manual <- function(seurat_obj, mapping) {
  Idents(seurat_obj) <- "seurat_clusters"
  clusters <- as.character(Idents(seurat_obj))
  celltype <- mapping[clusters]
  celltype[is.na(celltype)] <- clusters[is.na(celltype)]
  names(celltype) <- names(clusters)
  seurat_obj[["celltype"]] <- factor(celltype)
  Idents(seurat_obj) <- "celltype"
  return(seurat_obj)
}

seurat_list[["GSE72056"]] <- map_clusters_manual(seurat_list[["GSE72056"]], cluster2celltype_72056)
seurat_list[["GSE115978"]] <- map_clusters_manual(seurat_list[["GSE115978"]], cluster2celltype_115978)

# Save UMAP plots
for (dataset_name in names(seurat_list)) {
  p <- DimPlot(seurat_list[[dataset_name]], reduction = "umap", label = FALSE) +
    ggtitle(paste(dataset_name, "- Cell Types"))
  
  ggsave(
    filename = file.path(FIGURES_DIR, paste0(dataset_name, "_UMAP_CellTypes.pdf")),
    plot = p, width = 7, height = 6
  )
}
