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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_SeuratFinal.rds")

SeuratObjectBN1614 <- subset(SeuratObject, subset = sample_id %in% "BN1614") 

MetadataBN1614 <- SeuratObjectBN1614@meta.data

saveRDS(SeuratObjectBN1614,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_SeuratFinal_BN1614.rds")

write.table(MetadataBN1614,"Files/MetadataBN1614.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

SeuratObjectBN1578 <- subset(SeuratObject, subset = sample_id %in% "BN1578")

MetadataBN1578 <- SeuratObjectBN1578@meta.data

saveRDS(SeuratObjectBN1578,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_SeuratFinal_BN1578.rds")

write.table(MetadataBN1578,"Files/MetadataBN1578.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

cellbarcodesdf <- data.frame(cbind(rownames(SeuratObject@meta.data)[300000:301000],colnames(SeuratObject@assays$RNA@counts)[300000:301000]))

write.table(cellbarcodesdf, file = "Files/FINAL_cell_barcodes_df.txt", sep = "\t", quote = FALSE)

