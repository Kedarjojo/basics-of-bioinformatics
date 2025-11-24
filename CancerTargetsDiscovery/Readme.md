# Candidate Gene Analysis Across 10+ scRNA-seq Datasets

## Overview
This R script shows an example of automating the analysis of candidate genes across multiple single-cell RNA-seq (scRNA-seq) datasets. 
It is designed for scaling studies where you have 10 or more datasets and want to systematically compare gene expression across cell populations, clusters, or experimental conditions.

## Key Features
- Load multiple scRNA-seq datasets stored in a structured directory.  
- Preprocess each dataset individually: normalization, variable feature selection, etc.  
- Integrate datasets using Seurat’s **FindIntegrationAnchors** (supports batch integration to manage memory).  
- Merge integrated data into a single Seurat object for downstream analyses.  
- Automate candidate gene visualization across datasets and clusters.  
- Optionally generate pseudo-bulk averages per cluster for easier comparison.  

