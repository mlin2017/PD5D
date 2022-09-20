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
library(Nebulosa)
library(scCustomize)

dir.create("Figures/NebulosaPlots")

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_SeuratFinal.rds")

Idents(SeuratObject) <- "seurat_clusters"

Idents(SeuratObject) <- "CellSubtypes"

Top_Diff_LEGenes_df <- read.delim("Files/Top_50_DE_LEGenes.tsv", sep = "\t")

customPalette <- colorRampPalette(c("#FDF6B5","#F9C783","#4B1D91"))(200)

NebulosaFigureMaker <- function(SeurObj,Genes){
  for (i in Genes) {
    TempNebulosa <- Plot_Density_Custom(SeuratObject, i, custom_palette = customPalette, pt.size = 0.01) +
  facet_wrap(.~SeuratObject$case)
    ggsave(TempNebulosa, filename = paste("Figures/NebulosaPlots/",i,"_Nebulosa_SplitPlot.pdf",sep = ""), device = "pdf", width = 10, height = 4, units = "in")
  }}

#NebulosaFigureMaker(SeuratObject,Top_Diff_LEGenes_df$Top_Diff_LEGenes)

SeuratObject_UMAP_Clusters_nolegend <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.25, repel = TRUE) +
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 4),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"),
              legend.position = "none")

ggsave(SeuratObject_UMAP_Clusters_nolegend, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_UMAP_Clusters_nolegendALT.pdf",sep=""), device = "pdf", width = 10, height = 10, units = "in")

