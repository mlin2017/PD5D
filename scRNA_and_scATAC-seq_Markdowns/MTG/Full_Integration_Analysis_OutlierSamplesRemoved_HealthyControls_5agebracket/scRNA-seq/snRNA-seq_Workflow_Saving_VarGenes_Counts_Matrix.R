#Workflow up to elbow plot to determine cutoff, and saving intermediate seurat
#object as .rds file in /n/scratch3/users/j/jap0606/batch1to8

#set.seed(100)

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
library(stringr)
library(reshape2)
library(sciplot)

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls_5agebracket/FullIntegrationOSR_HealthyControls_5agebracket_MTG_PostClusterQCUMAP.rds")

VarSeuratObject <- SeuratObject[rownames(SeuratObject@assays$RNA) %in% SeuratObject@assays$RNA@var.features,]

VarSeuratCountsMatrix <- VarSeuratObject@assays$RNA@counts

saveRDS(VarSeuratCountsMatrix,"Files/HealthyControls_5AB_VarGenes_Count_Matrix.rds")

