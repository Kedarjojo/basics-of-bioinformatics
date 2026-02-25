source("scripts/config.R")

# Load raw file
df_72056 <- read.csv(
  "/home/joshi.ked/GSE72056.csv",
  check.names = FALSE
)

# Make rownames unique (not strictly necessary, but okay)
rownames(df_72056) <- make.unique(rownames(df_72056))

# Extract metadata (first 3 rows)
metadata_72056 <- df_72056[1:3, ]
rownames(metadata_72056) <- metadata_72056$Cell
metadata_72056 <- metadata_72056[, -1]

malignant_72056 <- as.factor(unlist(
  metadata_72056["malignant(1=no,2=yes,0=unresolved)", ]
))

celltype_72056 <- as.factor(unlist(
  metadata_72056["non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)", ]
))

tumor_72056 <- unlist(metadata_72056["tumor", ])

# Extract expression matrix (rows after metadata)
expr_df <- df_72056[-(1:3), ]

# Convert expression values to numeric
expr_numeric <- expr_df
expr_numeric[, -1] <- lapply(expr_numeric[, -1], as.numeric)

unique_genes <- make.unique(expr_numeric$Cell)
rownames(expr_numeric) <- unique_genes

# Drop the gene column (already in rownames)
expr_numeric <- expr_numeric[, -1]

# Convert to matrix
expr_matrix <- as.matrix(expr_numeric)

# Create metadata frame
meta72056 <- data.frame(
  Malignant = malignant_72056,
  CellType  = celltype_72056,
  Tumor     = tumor_72056,
  row.names = colnames(expr_matrix)
)

# Create Seurat object
seurat_72056 <- CreateSeuratObject(
  counts = expr_matrix,
  meta.data = meta72056
)

seurat_72056$dataset <- "GSE72056"

saveRDS(seurat_72056,
        file.path(OBJECTS_DIR, "seurat_72056.rds"))


expr_matrix_115978 <- read.csv(
  file.path(DATA_RAW, "GSE115978_counts.csv"),
  row.names = 1, check.names = FALSE
)

expr_matrix_115978 <- as.matrix(expr_matrix_115978)
mode(expr_matrix_115978) <- "numeric"

meta115978 <- read.csv(
  file.path(DATA_RAW, "GSE115978_cell.annotations.csv"),
  stringsAsFactors = FALSE
)

meta115978 <- meta115978[
  match(colnames(expr_matrix_115978), meta115978$cells), ]

rownames(meta115978) <- meta115978$cells
meta115978 <- meta115978[, -1]

seurat_115978 <- CreateSeuratObject(
  counts = expr_matrix_115978,
  meta.data = meta115978
)

seurat_115978$dataset <- "GSE115978"

saveRDS(seurat_115978,
        file.path(OBJECTS_DIR, "seurat_115978.rds"))
