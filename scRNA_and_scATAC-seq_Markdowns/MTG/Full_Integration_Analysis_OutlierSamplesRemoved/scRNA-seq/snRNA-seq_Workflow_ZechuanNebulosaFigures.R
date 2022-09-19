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
    ggsave(TempNebulosa, filename = paste("Figures/Zechuan_Figures/ALT",i,"_Nebulosa.pdf",sep = ""), device = "pdf", width = 4.8, height = 4, units = "in")
  }}

NebulosaFigureMaker(SeuratObject,MarkerGenes)






