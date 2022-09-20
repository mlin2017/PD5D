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
library(scales)
library(Nebulosa)
library(scCustomize)

dir.create("Figures/Zechuan_Figures")

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_SeuratFinal.rds")

Idents(SeuratObject) <- "seurat_clusters"

SeuratObject <- subset(SeuratObject, subset = MajorCellTypes != "Unknown_Cluster_66")

SeuratObject <- RunUMAP(SeuratObject, reduction = "harmony", dims = 1:60)

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons", `2` = "Oligodendrocytes", `3` = "GLU_Neurons", `4` = "GLU_Neurons", `5` = "GLU_Neurons", `6` = "Astrocytes", `7` = "GLU_Neurons", `8` = "GLU_Neurons",`9` = "GLU_Neurons",`10` = "GABA_Neurons", `11` = "GLU_Neurons", `12` = "Microglia",`13` = "GABA_Neurons",`14` = "OPCs",`15` = "GABA_Neurons", `16`= "GABA_Neurons", `17`="GABA_Neurons", `18`="Endothelial_Cells", `19`="Endothelial_Cells", `21` = "GABA_Neurons", `22` = "GLU_Neurons", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GLU_Neurons",`27` = "GLU_Neurons", `30` = "GABA_Neurons", `32` = "GLU_Neurons", `33` = "GLU_Neurons",`34` = "GABA_Neurons", `35` = "GLU_Neurons", `36` = "GABA_Neurons", `38` = "GABA_Neurons", `41` = "GABA_Neurons", `42` = "GABA_Neurons",`43` = "GABA_Neurons", `44` = "GLU_Neurons", `45` = "GLU_Neurons", `47` = "Unknown Cluster 8", `48` = "GABA_Neurons", `50` = "GABA_Neurons", `54` = "GLU_Neurons")

SeuratObject@meta.data$MajorCellTypes <- Idents(SeuratObject)

SeuratObject@meta.data$MajorCellTypes <- gsub("_"," ",SeuratObject@meta.data$MajorCellTypes)

Idents(SeuratObject) <- "MajorCellTypes"

levels(SeuratObject) <- c("GLU Neurons","GABA Neurons","Oligodendrocytes","Astrocytes","OPCs","Microglia","Endothelial Cells","Unknown Cluster 8")

MarkerGenes <- c("RBFOX3","SLC17A7","SLC17A6","GAD1","GAD2","PLP1","MBP","AQP4","GFAP","VCAN","BCAN","CX3CR1","P2RY12","FLT1","CLDN5")

ReducedMarkerGenes <- c("RBFOX3","SLC17A7","GAD2","PLP1","AQP4","VCAN","P2RY12","CLDN5")

##################################################################################################################################################################################

#UMAP Plots

FinalUMAPPlot_ZechuanFigures <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE, cols = c('GLU Neurons' = '#F8766D', 'Oligodendrocytes' = '#C49A00', 'GABA Neurons' = '#53B400', 'Astrocytes' = '#00C094', 'Microglia' = '#00B6EB', 'Endothelial Cells' = '#A58AFF', 'OPCs' = '#FB61D7', 'Unknown Cluster 8' = 'grey')) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(FinalUMAPPlot_ZechuanFigures, filename = "Figures/Zechuan_Figures/FinalUMAPPlot_Square.pdf", device = "pdf", width = 6, height = 6, units = "in")

ggsave(FinalUMAPPlot_ZechuanFigures, filename = "Figures/Zechuan_Figures/FinalUMAPPlot_Long.pdf", device = "pdf", width = 8, height = 6, units = "in")

FinalUMAPPlot_ZechuanFigures_nolabel <- DimPlot(SeuratObject, reduction = "umap", label = FALSE, pt.size = 0.01, label.size=2.5, repel = TRUE, cols = c('GLU Neurons' = '#F8766D', 'Oligodendrocytes' = '#C49A00', 'GABA Neurons' = '#53B400', 'Astrocytes' = '#00C094', 'Microglia' = '#00B6EB', 'Endothelial Cells' = '#A58AFF', 'OPCs' = '#FB61D7', 'Unknown Cluster 8' = 'grey')) +   
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(FinalUMAPPlot_ZechuanFigures_nolabel, filename = "Figures/Zechuan_Figures/FinalUMAPPlot_nolabel_Square.pdf", device = "pdf", width = 6, height = 6, units = "in")

ggsave(FinalUMAPPlot_ZechuanFigures_nolabel, filename = "Figures/Zechuan_Figures/FinalUMAPPlot_nolabel_Long.pdf", device = "pdf", width = 8, height = 6, units = "in")

FinalUMAPPlot_ZechuanFigures_nolegend <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE, cols = c('GLU Neurons' = '#F8766D', 'Oligodendrocytes' = '#C49A00', 'GABA Neurons' = '#53B400', 'Astrocytes' = '#00C094', 'Microglia' = '#00B6EB', 'Endothelial Cells' = '#A58AFF', 'OPCs' = '#FB61D7', 'Unknown Cluster 8' = 'grey')) +
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"),
        legend.position = "none")

ggsave(FinalUMAPPlot_ZechuanFigures_nolegend, filename = "Figures/Zechuan_Figures/FinalUMAPPlot_nolegend_Square.pdf", device = "pdf", width = 6, height = 6, units = "in")

ggsave(FinalUMAPPlot_ZechuanFigures_nolegend, filename = "Figures/Zechuan_Figures/FinalUMAPPlot_nolegend_Long.pdf", device = "pdf", width = 8, height = 6, units = "in")

########################################################################################################################################################################################

#Stacked Violin Plots

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(size = 12, angle = 90, face = "bold",hjust=0.95,vjust=0.2), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


MarkerViolinStack <- StackedVlnPlot(SeuratObject, features = MarkerGenes)

ggsave(MarkerViolinStack, filename = paste("Figures/Zechuan_Figures/MarkerViolinPlotStack.pdf",sep = ""), device = "pdf", width = 6, height = 16, units = "in")

ReducedMarkerViolinStack <- StackedVlnPlot(SeuratObject, features = ReducedMarkerGenes)

ggsave(ReducedMarkerViolinStack, filename = paste("Figures/Zechuan_Figures/MarkerSubsetViolinPlotStack.pdf",sep = ""), device = "pdf", width = 6, height = 10, units = "in")

#################################################################################################################################################################################################

#Nebulosa Plots

#Nebulosa_MarkerGenes <- plot_density(SeuratObject, features = MarkerGenes)

#ggsave(Nebulosa_MarkerGenes, filename = paste("Figures/Zechuan_Figures/Nebulosa_MarkerGenes.pdf",sep = ""), device = "pdf", width = 12, height = 10, units = "in")

#Nebulosa_ReducedMarkerGenes <- plot_density(SeuratObject, features = ReducedMarkerGenes)

#ggsave(Nebulosa_ReducedMarkerGenes, filename = paste("Figures/Zechuan_Figures/Nebulosa_ReducedMarkerGenes.pdf",sep = ""), device = "pdf", width = 10, height = 8, units = "in")

#customPalette <- colorRampPalette(c("#FDF6B5","#F9C783","#4B1D91"))(200)

customPalette <- colorRampPalette(c("#DEDDE6","#3C2692"))(200)

NebulosaFigureMaker <- function(SeurObj,Genes){
  for (i in Genes) {
    TempNebulosa <- Plot_Density_Custom(SeuratObject, i, custom_palette = customPalette, pt.size = 0.01)
    ggsave(TempNebulosa, filename = paste("Figures/Zechuan_Figures/",i,"_Nebulosa.pdf",sep = ""), device = "pdf", width = 4.8, height = 4, units = "in")
  }}

NebulosaFigureMaker(SeuratObject,MarkerGenes)






