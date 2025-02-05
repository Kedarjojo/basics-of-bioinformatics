### Create a main Seurat object with entire mouse brain (cereb/cortex/striat) data


genename<-read.delim("datasource/snRNA_MAPT_Mouse_Brain/aggregate_analysis_filtered/aggregate_analysis_cortex_filtered/filtered_feature_bc_matrix/features.tsv.gz", header = FALSE, sep = "\t") 
sample<-read.delim("datasource/snRNA_MAPT_Mouse_Brain/aggregate_analysis_filtered/aggregate_analysis_cortex_filtered/filtered_feature_bc_matrix/barcodes.tsv.gz", header = FALSE, sep = "\t") 
wt_data <- readMM(file = "datasource/snRNA_MAPT_Mouse_Brain/aggregate_analysis_filtered/aggregate_analysis_cortex_filtered/filtered_feature_bc_matrix/matrix.mtx.gz") 


rownames(wt_data)<-genename$V2 
colnames(wt_data)<-sample$V1 

meta1<-read.csv("datasource/snRNA_MAPT_Mouse_Brain/custom_analysis/mouse_cortex_snrna_01/integrated/csv_for_cellranger_and_loupe/barcode_vs_cell_type_for_loupe.csv",sep = ",") 
rownames(meta1)<-meta1$barcode 
meta2<-read.csv("datasource/snRNA_MAPT_Mouse_Brain/custom_analysis/mouse_cortex_snrna_01/integrated/csv_for_cellranger_and_loupe/barcode_vs_sample_name_for_loupe.csv",sep = ",") 
rownames(meta2)<-meta2$barcode
meta<-cbind(meta1,meta2) 

wt_data <- wt_data[!duplicated(rownames(wt_data)), ]

CNS_mice_brain <- CreateSeuratObject(counts = wt_data,meta.data = meta) 
CNS_mice_brain <- NormalizeData(CNS_mice_brain, normalization.method = "LogNormalize", scale.factor = 10000) 
CNS_mice_brain <- FindVariableFeatures(CNS_mice_brain, selection.method = "vst", nfeatures = 2000) 
all.genes <- rownames(CNS_mice_brain) 
CNS_mice_brain <- ScaleData(CNS_mice_brain, features = all.genes) 

CNS_mice_brain <- RunPCA(CNS_mice_brain, features = VariableFeatures(object = CNS_mice_brain)) 
CNS_mice_brain <- FindNeighbors(CNS_mice_brain, dims = 1:25) 
CNS_mice_brain <- FindClusters(CNS_mice_brain, resolution = 0.30)
CNS_mice_brain <- RunHarmony(CNS_mice_brain, "sample_name")
CNS_mice_brain <- RunUMAP(CNS_mice_brain, reduction = "harmony", dims =1:25)


##############################################################################################################################

## Attach the part of brain data you want to process to the main data

human_matrix <- ReadMtx(mtx = "datasource_kedar/MAPT_Human/Cereb10_11/outs/filtered_feature_bc_matrix/matrix.mtx", 
                        cells = "datasource_kedar/MAPT_Human/Cereb10_11/outs/filtered_feature_bc_matrix/barcodes.tsv",                        
                        features = "datasource_kedar/MAPT_Human/Cereb10_11/outs/filtered_feature_bc_matrix/features.tsv") 

mouse_matrix <- ReadMtx(mtx = "datasource_kedar/MAPT_Mouse/Cereb10_11/outs/filtered_feature_bc_matrix/matrix.mtx", 
                        cells = "datasource_kedar/MAPT_Mouse/Cereb10_11/outs/filtered_feature_bc_matrix/barcodes.tsv",                        
                        features = "datasource_kedar/MAPT_Mouse/Cereb10_11/outs/filtered_feature_bc_matrix/features.tsv") 


human_barcodes <- read.delim("datasource_kedar/MAPT_Human/Cereb10_11/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE) 
mouse_barcodes <- read.delim("datasource_kedar/MAPT_Mouse/Cereb10_11/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE) 



colnames(human_matrix) <- human_barcodes$V1 
colnames(mouse_matrix) <- mouse_barcodes$V1 

brain_barcodes <- colnames(CNS_mice_brain)

matched_human <- match(brain_barcodes, colnames(human_matrix)) 
matched_mouse <- match(brain_barcodes, colnames(mouse_matrix)) 

matched_human <- matched_human[!is.na(matched_human)] 
matched_mouse <- matched_mouse[!is.na(matched_mouse)] 

human_matrix <- human_matrix[, matched_human, drop = FALSE] 
mouse_matrix <- mouse_matrix[, matched_mouse, drop = FALSE] 

rownames(human_matrix) <- "human.Tfrc" 
rownames(mouse_matrix) <- "mouse.Tfrc" 

human_barcodes <- colnames(human_matrix) 
mouse_barcodes <- colnames(mouse_matrix) 

missing_human_barcodes <- setdiff(brain_barcodes, human_barcodes) 
missing_mouse_barcodes <- setdiff(brain_barcodes, mouse_barcodes) 

if (length(missing_human_barcodes) > 0) { 
  zero_human_matrix <- Matrix(0, nrow = 1, ncol = length(missing_human_barcodes), sparse = TRUE) 
  colnames(zero_human_matrix) <- missing_human_barcodes 
  rownames(zero_human_matrix) <- "human.Tfrc" 
  human_matrix <- cbind(human_matrix, zero_human_matrix) 
} 

if (length(missing_mouse_barcodes) > 0) { 
  zero_mouse_matrix <- Matrix(0, nrow = 1, ncol = length(missing_mouse_barcodes), sparse = TRUE) 
  colnames(zero_mouse_matrix) <- missing_mouse_barcodes 
  rownames(zero_mouse_matrix) <- "mouse.Tfrc" 
  mouse_matrix <- cbind(mouse_matrix, zero_mouse_matrix) 
} 

human_matrix <- human_matrix[, brain_barcodes, drop = FALSE] 
mouse_matrix <- mouse_matrix[, brain_barcodes, drop = FALSE] 

human_assay <- CreateAssayObject(counts = human_matrix) 
CNS_mice_brain_Cortex[["human.Tfrc"]] <- human_assay 

mouse_assay <- CreateAssayObject(counts = mouse_matrix) 

CNS_mice_brain_Cortex[["mouse.Tfrc"]] <- mouse_assay 


CNS_mice_brain_Cortex <- NormalizeData(CNS_mice_brain_Cortex, assay = "human.Tfrc") 
CNS_mice_brain_Cortex <- NormalizeData(CNS_mice_brain_Cortex,assay = "mouse.Tfrc")

CNS_mice_brain_Cortex <- ScaleData(CNS_mice_brain_Cortex, assay = "human.Tfrc") 
CNS_mice_brain_Cortex <- ScaleData(CNS_mice_brain_Cortex,assay = "mouse.Tfrc")

CNS_mice_brain_Cortex <- RunPCA(CNS_mice_brain_Cortex, features = VariableFeatures(object = CNS_mice_brain_Cortex))

CNS_mice_brain_Cortex <- FindNeighbors(CNS_mice_brain_Cortex, dims = 1:25) 

CNS_mice_brain_Cortex <- FindClusters(CNS_mice_brain_Cortex, resolution = 0.30)

CNS_mice_brain_Cortex <- RunHarmony(CNS_mice_brain_Cortex, "sample_name")

CNS_mice_brain_Cortex <- RunUMAP(CNS_mice_brain_Cortex, reduction = "harmony", dims =1:25)

DimPlot(CNS_mice_brain_Cortex, reduction = "umap", label = TRUE, repel = TRUE) 


Idents(CNS_mice_brain_Cortex)<-CNS_mice_brain_Cortex$predicted.per_sample.cell_type_detail_1

unique(CNS_mice_brain_Cortex$sample_name)

CNS_mice_brain_subset <- subset(CNS_mice_brain_Cortex, subset = sample_name == "Cortex_Veh")
