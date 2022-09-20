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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls/FullIntegrationOSR_HealthyControls_MTG_PostClusterQCUMAP.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons", `2` = "Oligodendrocytes", `3` = "GLU_Neurons", `4` = "GLU_Neurons", `5` = "GLU_Neurons", `6` = "Astrocytes", `7` = "GLU_Neurons", `8` = "GLU_Neurons",`9` = "GLU_Neurons",`10` = "GABA_Neurons", `11` = "Microglia", `12` = "GLU_Neurons",`13` = "GLU_Neurons",`14` = "GABA_Neurons",`15` = "Endothelial_Cells", `16`= "OPCs", `17`="GABA_Neurons", `18` = "GLU_Neurons", `19`="GABA_Neurons", `20`="Endothelial_Cells", `21` = "GABA_Neurons", `23` = "GLU_Neurons", `24` = "GABA_Neurons", `26` = "GABA_Neurons",`27` = "GLU_Neurons", `28` = "GLU_Neurons", `29` = "GLU_Neurons", `31` = "GLU_Neurons", `32` = "GABA_Neurons",`34` = "GABA_Neurons", `35` = "GABA_Neurons", `37` = "GABA_Neurons")

SeuratObject@meta.data$MajorCellTypes <- Idents(SeuratObject)

Idents(SeuratObject) <- "seurat_clusters"

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons_1", `2` = "Oligodendrocytes", `3` = "GLU_Neurons_2", `4` = "GLU_Neurons_3", `5` = "GLU_Neurons_4", `6` = "Astrocytes", `7` = "GLU_Neurons_5", `8` = "GLU_Neurons_6",`9` = "GLU_Neurons_7",`10` = "GABA_Neurons_1", `11` = "Microglia", `12` = "GLU_Neurons_8",`13` = "GLU_Neurons_9",`14` = "GABA_Neurons_2",`15` = "Endothelial_Cells_1", `16`= "OPCs", `17`="GABA_Neurons_3", `18` = "GLU_Neurons_10", `19`="GABA_Neurons_4", `20`="Endothelial_Cells_2", `21` = "GABA_Neurons_5", `23` = "GLU_Neurons_11", `24` = "GABA_Neurons_6", `26` = "GABA_Neurons_7",`27` = "GLU_Neurons_12", `28` = "GLU_Neurons_13", `29` = "GLU_Neurons_14", `31` = "GLU_Neurons_15", `32` = "GABA_Neurons_8",`34` = "GABA_Neurons_9", `35` = "GABA_Neurons_10", `37` = "GABA_Neurons_11")

SeuratObject@meta.data$CellSubtypes <- Idents(SeuratObject)

DietSeuratObject <- DietSeurat(SeuratObject, dimreducs = c("glmpca","harmony","umap"))

saveRDS(DietSeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls/FullIntegrationOSR_HC_MTG_DietSeuratFinal.rds")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls/FullIntegrationOSR_HC_MTG_SeuratFinal.rds")
