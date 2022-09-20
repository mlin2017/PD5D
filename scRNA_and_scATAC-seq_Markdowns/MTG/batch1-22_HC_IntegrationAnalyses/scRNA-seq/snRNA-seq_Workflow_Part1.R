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
library(glmpca)
library(SeuratWrappers)

dir.create("Files")
dir.create("Figures")

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/batch1-22_HC/Batch1-22_MTG_HC_CombinedFilteredMatrix.rds")

##########################################################################################################################

SeuratObject <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(SeuratObject)

SeuratObject <- ScaleData(SeuratObject, features = all.genes, verbose = FALSE)

#finding the top 30 principal components for cells
SeuratObject <- RunGLMPCA(SeuratObject, features=SeuratObject@assays$RNA@var.features, L = 50)

SeuratObject <- RunHarmony(SeuratObject, group.by.vars = c("sample_id","batch","sex","age_bracket"), plot_convergence = TRUE, reduction = "glmpca", theta = c(0.5,0.5,0.5,0.5))

harmony_embeddings <- Embeddings(SeuratObject, 'harmony')

SeuratObject_ElbowPlot <- ElbowPlot(SeuratObject, reduction = "harmony", ndims = 75)

ggsave2("Figures/SeuratObject_ElbowPlot.pdf", SeuratObject_ElbowPlot, device = "pdf", width = 4, height = 4, units = "in")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/batch1-22_HC/Batch1to22_MTG_HC_Part1.rds")


