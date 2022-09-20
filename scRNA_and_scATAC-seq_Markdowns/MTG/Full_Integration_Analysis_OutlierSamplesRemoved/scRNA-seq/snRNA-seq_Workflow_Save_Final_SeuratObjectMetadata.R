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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_SeuratFinal.rds")

metadata <- SeuratObject@meta.data

metadata$barcode <- colnames(SeuratObject@assays$RNA@counts)

write.table(metadata, file = "Files/FullIntegrationSeuratObject_FinaMetadataPlusBarcodes.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
