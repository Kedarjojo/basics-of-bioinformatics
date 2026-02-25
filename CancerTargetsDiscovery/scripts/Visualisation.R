source("scripts/config.R")

combined <- readRDS(
  file.path(OBJECTS_DIR, "combined_integrated.rds")
)

avg_expr <- AverageExpression(
  combined,
  assays = "RNA",
  features = CANDIDATE_GENES,
  group.by = "dataset"
)

png(file.path(FIGURES_DIR, "candidate_gene_heatmap.png"),
    width = 1000, height = 800)

pheatmap(avg_expr$RNA,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         main = "Candidate Gene Expression")

dev.off()