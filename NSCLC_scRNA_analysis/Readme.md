# NSCLC Single-Cell RNA-seq Analysis Pipeline

## Overview
This repository provides a complete **single-cell RNA-seq (scRNA-seq) analysis pipeline** for human **non-small cell lung cancer (NSCLC)** using tumor and matched normal lung tissue. The pipeline integrates multiple tools to explore cellular heterogeneity, identify cell types, detect marker genes, and perform trajectory analysis.

### Key Features of the Pipeline
- Data QC, normalization, and scaling with **Seurat**
- Doublet detection using **DoubletFinder**
- Cell type annotation using **SingleR** (human reference)
- Marker gene identification
- Trajectory inference using **Monocle3**
- Visualizations including **UMAPs**, **FeaturePlots**, and **ClusterPlots**

## Dataset
We use single-cell RNA-seq data from NSCLC patients (**LUAD** and **LUSC**) with both tumor and normal tissue samples. The data was generated using **10x Genomics Chromium V2 chemistry**.

### Selected Samples for Analysis
- **LUAD tumor-normal pairs:** 558T/N, 564T/N, 569T/N, 570T/N, 571T, 572T/N  
- **LUSC tumor-normal pairs:** 584TN1/2  
- **Additional LUAD tumor-normal pairs:** 593TN, 596TN, 630TN  

## Objective
The pipeline aims to:
1. Integrate multiple single-cell datasets to create a harmonized Seurat object
2. Perform QC and doublet removal
3. Identify and annotate major cell types
4. Detect marker genes for clusters and cell types
5. Explore cellular trajectories and differentiation using Monocle3

## Citation
Leader AM, Grout JA, Maier BB, Nabet BY, et al.  
*Single-cell analysis of human non-small cell lung cancer lesions refines tumor classification and patient stratification.*  
**Cancer Cell** 2021 Dec 13;39(12):1594-1609.e12.  
PMID: [34767762](https://pubmed.ncbi.nlm.nih.gov/34767762/)
