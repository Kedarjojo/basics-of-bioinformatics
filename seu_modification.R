human_matrix <- ReadMtx(mtx = "~/Desktop/Data/CNS_Delivery_workstream_6/snRNA-seq_Cortex_Mice/Data/scRNA-Mouse_Brain_CNS_W6/TFRC_Human_and_Mouse/Human/Veh/outs/filtered_feature_bc_matrix/matrix.mtx.gz", 
                        cells = "~/Desktop/Data/CNS_Delivery_workstream_6/snRNA-seq_Cortex_Mice/Data/scRNA-Mouse_Brain_CNS_W6/TFRC_Human_and_Mouse/Human/Veh/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",                        
                        features = "~/Desktop/Data/CNS_Delivery_workstream_6/snRNA-seq_Cortex_Mice/Data/scRNA-Mouse_Brain_CNS_W6/TFRC_Human_and_Mouse/Human/Veh/outs/filtered_feature_bc_matrix/features.tsv.gz") 

mouse_matrix <- ReadMtx(mtx = "~/Desktop/Data/CNS_Delivery_workstream_6/snRNA-seq_Cortex_Mice/Data/scRNA-Mouse_Brain_CNS_W6/TFRC_Human_and_Mouse/Mouse/Veh/outs/filtered_feature_bc_matrix/matrix.mtx.gz",                        
                        cells = "~/Desktop/Data/CNS_Delivery_workstream_6/snRNA-seq_Cortex_Mice/Data/scRNA-Mouse_Brain_CNS_W6/TFRC_Human_and_Mouse/Mouse/Veh/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                        features = "~/Desktop/Data/CNS_Delivery_workstream_6/snRNA-seq_Cortex_Mice/Data/scRNA-Mouse_Brain_CNS_W6/TFRC_Human_and_Mouse/Mouse/Veh/outs/filtered_feature_bc_matrix/features.tsv.gz") 

human_barcodes <- read.delim("~/Desktop/Data/CNS_Delivery_workstream_6/snRNA-seq_Cortex_Mice/Data/scRNA-Mouse_Brain_CNS_W6/TFRC_Human_and_Mouse/Human/Veh/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", header = FALSE) 

mouse_barcodes <- read.delim("~/Desktop/Data/CNS_Delivery_workstream_6/snRNA-seq_Cortex_Mice/Data/scRNA-Mouse_Brain_CNS_W6/TFRC_Human_and_Mouse/Mouse/Veh/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", header = FALSE) 

colnames(human_matrix) <- human_barcodes$V1 

colnames(mouse_matrix) <- mouse_barcodes$V1 

existing_barcodes <- colnames(subset_veh) 

matched_human <- match(existing_barcodes, colnames(human_matrix)) 
matched_mouse <- match(existing_barcodes, colnames(mouse_matrix)) 

matched_human <- matched_human[!is.na(matched_human)] matched_mouse <- matched_mouse[!is.na(matched_mouse)] human_matrix <- human_matrix[, matched_human, drop = FALSE] 
mouse_matrix <- mouse_matrix[, matched_mouse, drop = FALSE] 
rownames(human_matrix) <- "human.Tfrc" 
rownames(mouse_matrix) <- "mouse.Tfrc" 
human_barcodes <- colnames(human_matrix) 
mouse_barcodes <- colnames(mouse_matrix) 
missing_human_barcodes <- setdiff(existing_barcodes, human_barcodes) 
missing_mouse_barcodes <- setdiff(existing_barcodes, mouse_barcodes) 
if (length(missing_human_barcodes) > 0) { 
  zero_human_matrix <- Matrix(0, nrow = 1, ncol = length(missing_human_barcodes), sparse = TRUE) 
  colnames(zero_human_matrix) <- missing_human_barcodes 
  rownames(zero_human_matrix) <- "human.Tfrc" 
  human_matrix <- cbind(human_matrix, zero_human_matrix) 
} 
if (length(missing_mouse_barcodes) > 0) { 
  zero_mouse_matrix <- Matrix(0, nrow = 1, ncol = length(missing_mouse_barcodes), sparse = TRUE) 
  colnames(zero_mouse_matrix) <- missing_mouse_barcodes 
  rownames(zero_mouse_matrix) <- "mouse.Tfrc" mouse_matrix <- cbind(mouse_matrix, zero_mouse_matrix) 
} 

human_matrix <- human_matrix[, existing_barcodes, drop = FALSE] 
mouse_matrix <- mouse_matrix[, existing_barcodes, drop = FALSE] 
human_assay <- CreateAssayObject(counts = human_matrix) 
subset_veh[["human.Tfrc"]] <- human_assay 

mouse_assay <- CreateAssayObject(counts = mouse_matrix) 
subset_veh[["mouse.Tfrc"]] <- mouse_assay 
subset_veh <- NormalizeData(subset_veh, assay = "human.Tfrc") 
subset_veh <- NormalizeData(subset_veh, assay = "mouse.Tfrc") 
subset_veh <- ScaleData(subset_veh, assay = "human.Tfrc") 
subset_veh <- ScaleData(subset_veh, assay = "mouse.Tfrc")