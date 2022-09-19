library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(dplyr)
library(psych)
library(pheatmap)
library(harmony)
#library(clusterProfiler)
#library(org.Hs.eg.db)
library(DOSE)
library(GOSemSim)
library(enrichplot)
library(glmpca)
library(SeuratWrappers)
library(stringr)

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls_5agebracket/FullIntegrationOSR_HealthyControls_5agebracket_MTG_PostClusterQCUMAP.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "Oligodendrocytes", `2` = "GLU_Neurons", `3` = "GLU_Neurons", `4` = "GLU_Neurons", `5` = "GLU_Neurons", `6` = "GLU_Neurons", `7` = "Astrocytes", `8` = "GLU_Neurons",`9` = "GABA_Neurons",`10` = "GABA_Neurons", `11` = "GLU_Neurons", `12` = "Microglia",`13` = "GLU_Neurons",`14` = "GABA_Neurons",`15` = "GLU_Neurons", `16`= "GLU_Neurons", `17`="OPCs", `19`="Endothelial_Cells", `20`="Endothelial_Cells", `21` = "GABA_Neurons", `22` = "GABA_Neurons", `23` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GABA_Neurons",`27` = "GLU_Neurons", `28` = "GLU_Neurons", `29` = "GLU_Neurons", `32` = "GABA_Neurons", `33` = "GLU_Neurons",`34` = "GLU_Neurons", `35` = "GLU_Neurons", `36` = "GABA_Neurons", `37` = "GABA_Neurons")

SeuratObject@meta.data$MajorCellTypes <- Idents(SeuratObject)

Idents(SeuratObject) <- "seurat_clusters"

SeuratObject <- RenameIdents(SeuratObject, `1` = "Oligodendrocytes", `2` = "GLU_Neurons_1", `3` = "GLU_Neurons_2", `4` = "GLU_Neurons_3", `5` = "GLU_Neurons_4", `6` = "GLU_Neurons_5", `7` = "Astrocytes", `8` = "GLU_Neurons_6",`9` = "GABA_Neurons_1",`10` = "GABA_Neurons_2", `11` = "GLU_Neurons_7", `12` = "Microglia",`13` = "GLU_Neurons_8",`14` = "GABA_Neurons_3",`15` = "GLU_Neurons_9", `16`= "GLU_Neurons_10", `17`="OPCs", `19`="Endothelial_Cells_1", `20`="Endothelial_Cells_2", `21` = "GABA_Neurons_4", `22` = "GABA_Neurons_5", `23` = "GLU_Neurons_11", `25` = "GLU_Neurons_12", `26` = "GABA_Neurons_6",`27` = "GLU_Neurons_13", `28` = "GLU_Neurons_14", `29` = "GLU_Neurons_15", `32` = "GABA_Neurons_7", `33` = "GLU_Neurons_16",`34` = "GLU_Neurons_17", `35` = "GLU_Neurons_18", `36` = "GABA_Neurons_8", `37` = "GABA_Neurons_9")

SeuratObject@meta.data$CellSubtypes <- Idents(SeuratObject)

DietSeuratObject <- DietSeurat(SeuratObject, dimreducs = c("glmpca","harmony","umap"))

saveRDS(DietSeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls_5agebracket/FullIntegrationOSR_HC_5agebracket_DietSeuratFinal.rds")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls_5agebracket/FullIntegrationOSR_HC_5agebracket_SeuratFinal.rds")

