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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_PostClusterQCUMAP.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons", `2` = "Oligodendrocytes", `3` = "GLU_Neurons", `4` = "GLU_Neurons", `5` = "GLU_Neurons", `6` = "Astrocytes", `7` = "GLU_Neurons", `8` = "GLU_Neurons",`9` = "GLU_Neurons",`10` = "GABA_Neurons", `11` = "GLU_Neurons", `12` = "Microglia",`13` = "GABA_Neurons",`14` = "OPCs",`15` = "GABA_Neurons", `16`= "GABA_Neurons", `17`="GABA_Neurons", `18`="Endothelial_Cells", `19`="Endothelial_Cells", `21` = "GABA_Neurons", `22` = "GLU_Neurons", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GLU_Neurons",`27` = "GLU_Neurons", `30` = "GABA_Neurons", `32` = "GLU_Neurons", `33` = "GLU_Neurons",`34` = "GABA_Neurons", `35` = "GLU_Neurons", `36` = "GABA_Neurons", `38` = "GABA_Neurons", `41` = "GABA_Neurons", `42` = "GABA_Neurons",`43` = "GABA_Neurons", `44` = "GLU_Neurons", `45` = "GLU_Neurons", `47` = "TEMRA_T_Cells", `48` = "GABA_Neurons", `50` = "GABA_Neurons", `54` = "GLU_Neurons", `66` = "Unknown_Cluster_66")

SeuratObject@meta.data$MajorCellTypes <- Idents(SeuratObject)

Idents(SeuratObject) <- "seurat_clusters"

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons_1", `2` = "Oligodendrocytes", `3` = "GLU_Neurons_2", `4` = "GLU_Neurons_3", `5` = "GLU_Neurons_4", `6` = "Astrocytes", `7` = "GLU_Neurons_5", `8` = "GLU_Neurons_6",`9` = "GLU_Neurons_7",`10` = "GABA_Neurons_1", `11` = "GLU_Neurons_8", `12` = "Microglia",`13` = "GABA_Neurons_2",`14` = "OPCs",`15` = "GABA_Neurons_3", `16`= "GABA_Neurons_4", `17`="GABA_Neurons_5", `18`="Endothelial_Cells_1", `19`="Endothelial_Cells_2", `21` = "GABA_Neurons_6", `22` = "GLU_Neurons_9", `23` = "GLU_Neurons_10", `24` = "GLU_Neurons_11", `25` = "GLU_Neurons_12", `26` = "GLU_Neurons_13",`27` = "GLU_Neurons_14", `30` = "GABA_Neurons_7", `32` = "GLU_Neurons_15", `33` = "GLU_Neurons_16",`34` = "GABA_Neurons_8", `35` = "GLU_Neurons_17", `36` = "GABA_Neurons_9", `38` = "GABA_Neurons_10", `41` = "GABA_Neurons_11", `42` = "GABA_Neurons_12",`43` = "GABA_Neurons_13", `44` = "GLU_Neurons_18", `45` = "GLU_Neurons_19", `47` = "TEMRA_T_Cells", `48` = "GABA_Neurons_14", `50` = "GABA_Neurons_15", `54` = "GLU_Neurons_20", `66` = "Unknown_Cluster_66")

SeuratObject@meta.data$CellSubtypes <- Idents(SeuratObject)

Idents(SeuratObject) <- "seurat_clusters"

DietSeuratObject <- DietSeurat(SeuratObject, dimreducs = c("glmpca","harmony","umap"))

metadata <- SeuratObject@meta.data

saveRDS(DietSeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_DietSeuratFinal.rds")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_SeuratFinal.rds")

write.table(metadata, file = "Files/FullIntegration_OutlierSamplesRemoved_FinaMetadata.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
