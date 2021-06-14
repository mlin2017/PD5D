#Workflow up to elbow plot to determine cutoff, and saving intermediate seurat
#object as .rds file in /n/scratch3/users/j/jap0606/batch1to8

library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(dplyr)
library(psych)
library(pheatmap)
library(harmony)
#library(clusterProfiler)
library(DOSE)
library(GOSemSim)
library(enrichplot)
library(stringr)

Batch1to8_MTG <- readRDS("Files/Batch1to8_MTG_preclustering.rds")

all.genes <- rownames(Batch1to8_MTG)
Batch1to8_MTG <- ScaleData(Batch1to8_MTG, features = all.genes, verbose = FALSE)

#finding the top 30 principal components for cells
Batch1to8_MTG <- RunPCA(Batch1to8_MTG, npcs = 30, verbose = FALSE)

Batch1to8_MTG <- RunHarmony(Batch1to8_MTG, group.by.vars = c("case","batch"), plot_convergence = TRUE)

harmony_embeddings <- Embeddings(Batch1to8_MTG, 'harmony')

Batch1to8_MTG_ElbowPlot <- ElbowPlot(Batch1to8_MTG)

ggsave2("Figures/Batch1to8_MTG_ElbowPlot.pdf", Batch1to8_MTG_ElbowPlot, device = "pdf", width = 4, height = 4, units = "in")

saveRDS(Batch1to8_MTG,"/n/scratch3/users/j/jap0606/batch1to8/Batch1to8_MTG_Part1.rds")
