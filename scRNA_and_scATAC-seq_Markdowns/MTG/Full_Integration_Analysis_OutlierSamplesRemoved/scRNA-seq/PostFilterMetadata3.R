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
library(tidyr)
library(scDblFinder)

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration/FullIntegrationOSR_MTG_PostClusterQCUMAP.rds")

metadataDF <- SeuratObject@meta.data

write.table(metadataDF, file = "Files/PostFilteringMetadata_Table3.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
