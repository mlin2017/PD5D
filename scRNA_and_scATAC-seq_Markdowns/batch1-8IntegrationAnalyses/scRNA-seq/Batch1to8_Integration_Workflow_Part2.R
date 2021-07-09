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
#library(org.Hs.eg.db)
library(DOSE)
library(GOSemSim)
library(enrichplot)
library(stringr)
library(reshape2)
library(sciplot)

Batch1to8_MTG <- readRDS("/n/scratch3/users/j/jap0606/batch1to8/Batch1to8_MTG_Part1.rds")

Batch1to8_MTG <- FindNeighbors(Batch1to8_MTG, reduction = "harmony", dims = 1:20)
Batch1to8_MTG <- FindClusters(Batch1to8_MTG, resolution = 0.5, algorithm = 4, method = "igraph")
Batch1to8_MTG <- RunUMAP(Batch1to8_MTG, reduction = "harmony", dims = 1:20)

saveRDS(Batch1to8_MTG,"/n/scratch3/users/j/jap0606/batch1to8/Batch1to8_MTG_Part2.rds")

Batch1to8_MTG_UMAP_Clusters <- DimPlot(Batch1to8_MTG, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(Batch1to8_MTG_UMAP_Clusters, filename = "Figures/Batch1to8_MTG_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 4, units = "in")

Batch1to8_MTG_Case_Group_UMAP_Clusters <- DimPlot(Batch1to8_MTG, reduction = "umap", group.by = "case", pt.size = 0.1, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(Batch1to8_MTG_Case_Group_UMAP_Clusters, filename = "Figures/Batch1to8_MTG_Case_Group_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 4, units = "in")

Batch1to8_MTG_Region_Group_UMAP_Clusters <- DimPlot(Batch1to8_MTG, reduction = "umap", group.by = "batch",pt.size = 0.1, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(Batch1to8_MTG_Region_Group_UMAP_Clusters, filename = "Figures/Batch1to8_MTG_Region_Group_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 4, units = "in")

Batch1to8_MTG_Case_Split_UMAP_Clusters <- DimPlot(Batch1to8_MTG, reduction = "umap", split.by = "case", label = TRUE, ncol = 1, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(Batch1to8_MTG_Case_Split_UMAP_Clusters, filename = "Figures/Batch1to8_MTG_Case_Split_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 8, units = "in")

Batch1to8_MTG_Region_Split_UMAP_Clusters <- DimPlot(Batch1to8_MTG, reduction = "umap", split.by = "batch", label = TRUE, ncol = 1, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(Batch1to8_MTG_Region_Split_UMAP_Clusters, filename = "Figures/Batch1to8_MTG_Region_Split_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 32, units = "in")

MarkerGenes <- c("ENO2","RBFOX3","SLC17A6","SLC17A7","SLC32A1","GAD1","GAD2","AQP4","GFAP","PLP1","MBP","VCAN","BCAN","CX3CR1","P2RY12","FLT1","CLDN5","IL7R","CD96","CD8A","RELN","CALB2","CNR1")

data_barplot <- FetchData(Batch1to8_MTG, vars = c("ident",MarkerGenes), slot = "data")

#data_barplot <- FetchData(Batch1to8_MTG, vars = c("ident",rownames(Batch1to8_MTG@assays$RNA@counts)), slot = "counts")

#scaledrowSums <- 1e6/rowSums(data_barplot[2:length(colnames(data_barplot))])

#swpt_barplot <- sweep(data_barplot[,2:length(colnames(data_barplot))],1,scaledrowSums,FUN = "*")

#swpt_data_barplot_markers <- swpt_barplot[,which(colnames(swpt_barplot) %in% MarkerGenes)]

#swpt_data_barplot_markers$ident <- as.vector(data_barplot$ident)

#data_barplot_melt <- melt(swpt_data_barplot_markers)

data_barplot_melt <- melt(data_barplot)

data_barplot_melt$ident <- as.vector(data_barplot_melt$ident)
data_barplot_melt$variable <- as.vector(data_barplot_melt$variable)
data_barplot_melt$value <- as.numeric(as.vector(data_barplot_melt$value))

data_barplot_melt_sum <- group_by(data_barplot_melt,ident,variable) %>% summarise(mean = mean(value), SE = se(value))

data_barplot_melt_sum$ident <- factor(data_barplot_melt_sum$ident, levels = unique(data_barplot_melt_sum$ident))

data_barplot_melt_sum$variable <- factor(data_barplot_melt_sum$variable, levels = unique(MarkerGenes))

Batch1to8_barchart <- ggplot(data_barplot_melt_sum, aes(x = ident, y = mean, fill = ident)) + 
        geom_bar(aes(x = ident, y = mean), stat = "identity", alpha = 1) + 
        geom_errorbar(aes(x = ident, ymin = mean-SE, ymax = mean+SE, colour = ident), width = 0.4, alpha = 0.9, size = 0.5) + 
        ggplot2::facet_grid(rows = vars(variable), scales = "free_y", switch = "y") + 
        theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 45, face = "bold", vjust = 0.5),
              axis.text.y = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
              strip.background = element_blank(), strip.placement = "outside", 
              strip.text.y = element_text(size = 12, angle = 180, face = "bold"),
              strip.text.y.left = element_text(angle = 0)) + NoLegend()


ggsave(Batch1to8_barchart,filename = "Files/Batch1to8_Marker_Barchart", device = "pdf", width = 12, height = 14, units = "in")


Markerggplots <- function(SeurObj,Genes){
  for (i in Genes) {
    TempViolin <- VlnPlot(SeurObj, features = i ,pt.size = 0)
    ggsave(TempViolin, filename = paste("Figures/",i,"_VlnPlot.pdf",sep = ""), device = "pdf", width = 12, height = 4, units = "in")
}}

Markerggplots(Batch1to8_MTG,MarkerGenes)



Markerggplotspt1 <- function(SeurObj,Genes){
  for (i in Genes) {
    TempViolin <- VlnPlot(SeurObj, features = i ,pt.size = 1)
    ggsave(TempViolin, filename = paste("Figures/",i,"_Pt1_VlnPlot.pdf",sep = ""), device = "pdf", width = 12, height = 4, units = "in")
}}

Markerggplotspt1(Batch1to8_MTG,MarkerGenes)


