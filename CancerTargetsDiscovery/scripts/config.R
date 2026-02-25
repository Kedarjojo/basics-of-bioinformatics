# ==============================
# Global Configuration
# ==============================

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)

# Project paths (relative)
BASE_DIR <- "."
DATA_RAW <- file.path(BASE_DIR, "data/raw")
DATA_PROCESSED <- file.path(BASE_DIR, "data/processed")
RESULTS_DIR <- file.path(BASE_DIR, "results")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR <- file.path(RESULTS_DIR, "tables")
OBJECTS_DIR <- file.path(RESULTS_DIR, "objects")

# Create directories if they don't exist
dir.create(DATA_PROCESSED, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Candidate genes
CANDIDATE_GENES <- c(
  "CASKIN2","EMC9","PDIK1L","DBNDD2",
  "FAM171A2","C1orf174","LOC124903857",
  "TMEM161B","ZNF808"
)