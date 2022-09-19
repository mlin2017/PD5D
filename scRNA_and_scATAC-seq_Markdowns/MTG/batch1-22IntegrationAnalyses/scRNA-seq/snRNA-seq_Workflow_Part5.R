#Workflow up to elbow plot to determine cutoff, and saving intermediate seurat
#object as .rds file in /n/scratch3/users/j/jap0606/batch1to8

set.seed(100)

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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/batch1-22/Batch1-22_MTG_Part2.rds")

DimPlot(SeuratObject, label = TRUE, repel = TRUE, pt.size = 0, label.size = 2.5)

ClusterExpMed <- function(ClusterIdent, SeuObj) {
  IdentObj <- subset(SeuObj, idents = ClusterIdent)
  MedExp <- apply(IdentObj@assays$RNA@data, 1, median)
  return(MedExp)
}

ExpVectorPerCluster <- lapply(sort(as.numeric(unique(SeuratObject@active.ident))),FUN = ClusterExpMed, SeuObj = SeuratObject)

ExpVectorPerClusterDF <- as.data.frame(do.call(cbind, ExpVectorPerCluster))

colnames(ExpVectorPerClusterDF) <- sort(unique(SeuratObject@active.ident))

PanglaoDB <- read.delim("Files/PanglaoDB_markers_27_Mar_2020.tsv", stringsAsFactors = FALSE)

PanglaoDBBrain <- PanglaoDB[PanglaoDB$organ %in% "Brain",]

PanglaoDBBrain <- PanglaoDBBrain[!PanglaoDBBrain$species %in% "Mm",]

#1 + sqrt((length(unique(PanglaoDBBrain$cell.type)) - sum(PanglaoDBBrain$official.gene.symbol == "TH"))/(length(unique(PanglaoDBBrain$cell.type)) - 1))

AssignGeneWeights <- function(Gene) {
  Weight <- 1 + sqrt((length(unique(PanglaoDBBrain$cell.type)) - sum(PanglaoDBBrain$official.gene.symbol %in% Gene))/(length(unique(PanglaoDBBrain$cell.type)) - 1))
  return(Weight)
}

GeneWeights <- unlist(lapply(unique(PanglaoDBBrain$official.gene.symbol), AssignGeneWeights))

GeneWeightsDF <- as.data.frame(cbind(unique(PanglaoDBBrain$official.gene.symbol), GeneWeights))

colnames(GeneWeightsDF) <- c("Gene","Weight")

PanglaoDBBrain <- PanglaoDBBrain[PanglaoDBBrain$official.gene.symbol %in% rownames(SeuratObject@assays$RNA@counts),]

GeneWeightsDF <- GeneWeightsDF[GeneWeightsDF$Gene %in% PanglaoDBBrain$official.gene.symbol,]

PanglaoDBBrain$GeneWeights <- GeneWeightsDF$Weight[match(PanglaoDBBrain$official.gene.symbol,GeneWeightsDF$Gene)]

PanglaoDBBrain[PanglaoDBBrain$official.gene.symbol == "TH",]

CellTypes <- unique(PanglaoDBBrain$cell.type)

PanglaoDBBrain$GeneWeights <- as.numeric(PanglaoDBBrain$GeneWeights)

CalculateCTA <- function(Cell, ClusterFrame, MarkerFrame) {
  CellMarkerFrame <- MarkerFrame[MarkerFrame$cell.type %in% Cell,]
  ClusterFrame$Genes <- rownames(ClusterFrame)
  CellClusterFrame <- as.data.frame(ClusterFrame[ClusterFrame$Genes %in% CellMarkerFrame$official.gene.symbol,])
  CellMarkerFrame <- CellMarkerFrame[order(match(CellMarkerFrame$official.gene.symbol,rownames(CellClusterFrame))),]
  CTAUnscaled <- CellClusterFrame[[1]]*CellMarkerFrame$GeneWeights
  CTAScore <- sum(CTAUnscaled/((length(CellMarkerFrame$official.gene.symbol))^1/3))
  return(CTAScore)
}

CellClusterAssignment=NULL
for (i in colnames(ExpVectorPerClusterDF)){
  ClusterColumn <- as.data.frame(ExpVectorPerClusterDF[[i]])
  rownames(ClusterColumn) <- rownames(ExpVectorPerClusterDF)
  colnames(ClusterColumn) <- "MedianExpression"
  CTAVector <- unlist(lapply(CellTypes, FUN = CalculateCTA, ClusterFrame = ClusterColumn, MarkerFrame = PanglaoDBBrain))
  names(CTAVector) <- CellTypes
  Assignment <- names(CTAVector[CTAVector == max(CTAVector)])
  CellClusterAssignment <- c(CellClusterAssignment,Assignment)
}

write.table(as.data.frame(CellClusterAssignment), file = "Files/CellClusterAssignment.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


SeuratObject <- RenameIdents(SeuratObject, `1` = "Pyramidal Cells", `2` = "Pyramidal Cells", `3` = "Oligodendrocytes", `4` = "Pyramidal Cells", `5` = "Glu Neurons", `6` = "Pyramidal Cells", `7` = "Pyramidal Cells", `8` = "Purkinje Cells",`9` = "Purkinje Cells",`10` = "Adrenergic Neurons", `11` = "Radial Glia Cells",`12` = "Pyramidal Cells",`13` = "Glu Neurons",`14` = "OPCs",`15` = "Pyramidal Cells", `16`= "Pyramidal Cells", `17`="Pyramidal Cells", `18`="GABA Neurons", `19`="Purkinje Cells",  `20`="Neural Stem Cells", `21` = "Satellite Glia Cells", `22` = "Glu Neurons", `23` = "Pyramidal Cells", `24` = "Pyramidal Cells", `25` = "Glu Neurons", `26` = "Pyramidal Cells",`27` = "Purkinje Cells",`28` = "Pyramidal Cells",`29` = "OPCs",`30` = "Glu Neurons")


UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(UMAP_Clusters, filename = paste("Figures/PanglaoDB_Auto_Assigned_",SeuratObject@project.name,"_MTG_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 4, units = "in")



