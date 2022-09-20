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
library(sciplot)
library(reshape2)
library(MAST)
library(fgsea)
library(SeuratWrappers)
library(purrr)

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls/FullIntegrationOSR_HC_MTG_SeuratFinal.rds")

Idents(SeuratObject) <- "MajorCellTypes"

MarkerParkinsonsGenes <- c("SNCA", "PRKN", "UCHL1", "PINK1", "PARK7", "LRRK2", "ATP13A2", "GIGYF2", "HTRA2", "PLA2G6", "FBXO7", "VPS35", "EIF4G1", "DNAJC6", "SYNJ1", "DNAJC13", "CHCHD2", "VPS13C", "GBA")

data_barplot <- FetchData(SeuratObject, vars = c("ident",MarkerParkinsonsGenes), slot = "data")
data_barplot[,2:ncol(data_barplot)] <- apply(as.matrix(data_barplot[,2:ncol(data_barplot)]),2,expm1)
data_barplot_melt <- melt(data_barplot)

data_barplot_melt$ident <- as.vector(data_barplot_melt$ident)
data_barplot_melt$variable <- as.vector(data_barplot_melt$variable)
data_barplot_melt$value <- as.numeric(as.vector(data_barplot_melt$value))

data_barplot_melt_sum <- group_by(data_barplot_melt,ident,variable) %>% summarise(mean = mean(value), SE = se(value))

data_barplot_melt_sum$ident <- factor(data_barplot_melt_sum$ident, levels = unique(data_barplot_melt_sum$ident))

data_barplot_melt_sum$variable <- factor(data_barplot_melt_sum$variable, levels = unique(MarkerParkinsonsGenes))

write.table(data_barplot_melt_sum, file = paste("Files/MajorCellTypes_Parkisons_Mendelian_Genes_BarChartTable.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

MonikaGenes <- c("SNCA", "ADRB1", "ADRB2", "ADRB3", "PM20D1", "ATP5F1A", "LARS2", "UCP2", "UCP3", "SLC25A27")

data_barplot <- FetchData(SeuratObject, vars = c("ident",MonikaGenes), slot = "data")
data_barplot[,2:ncol(data_barplot)] <- apply(as.matrix(data_barplot[,2:ncol(data_barplot)]),2,expm1)
data_barplot_melt <- melt(data_barplot)

data_barplot_melt$ident <- as.vector(data_barplot_melt$ident)
data_barplot_melt$variable <- as.vector(data_barplot_melt$variable)
data_barplot_melt$value <- as.numeric(as.vector(data_barplot_melt$value))

data_barplot_melt_sum <- group_by(data_barplot_melt,ident,variable) %>% summarise(mean = mean(value), SE = se(value))

data_barplot_melt_sum$ident <- factor(data_barplot_melt_sum$ident, levels = unique(data_barplot_melt_sum$ident))

data_barplot_melt_sum$variable <- factor(data_barplot_melt_sum$variable, levels = unique(MonikaGenes))

write.table(data_barplot_melt_sum, file = paste("Files/MajorCellTypes_Monika_Genes_BarChartTable.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

BeaGenes <- c("RIMS2", "RIMS1", "SNCA", "BIN3", "BIN1", "SH3GL2")

data_barplot <- FetchData(SeuratObject, vars = c("ident",BeaGenes), slot = "data")
data_barplot[,2:ncol(data_barplot)] <- apply(as.matrix(data_barplot[,2:ncol(data_barplot)]),2,expm1)
data_barplot_melt <- melt(data_barplot)

data_barplot_melt$ident <- as.vector(data_barplot_melt$ident)
data_barplot_melt$variable <- as.vector(data_barplot_melt$variable)
data_barplot_melt$value <- as.numeric(as.vector(data_barplot_melt$value))

data_barplot_melt_sum <- group_by(data_barplot_melt,ident,variable) %>% summarise(mean = mean(value), SE = se(value))

data_barplot_melt_sum$ident <- factor(data_barplot_melt_sum$ident, levels = unique(data_barplot_melt_sum$ident))

data_barplot_melt_sum$variable <- factor(data_barplot_melt_sum$variable, levels = unique(BeaGenes))

write.table(data_barplot_melt_sum, file = paste("Files/MajorCellTypes_Bea_Genes_BarChartTable.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

#BeaGenes <- c("RIMS2", "RIMS1", "SNCA", "BIN3", "BIN1", "SH3GL2")

#data_barplot <- FetchData(SeuratObject, vars = c("ident",BeaGenes), slot = "data")
#data_barplot[,2:ncol(data_barplot)] <- apply(as.matrix(data_barplot[,2:ncol(data_barplot)]),2,expm1)
#data_barplot_melt <- melt(data_barplot)

#data_barplot_melt$ident <- as.vector(data_barplot_melt$ident)
#data_barplot_melt$variable <- as.vector(data_barplot_melt$variable)
#data_barplot_melt$value <- as.numeric(as.vector(data_barplot_melt$value))

#data_barplot_melt_sum <- group_by(data_barplot_melt,ident,variable) %>% summarise(mean = mean(value), SE = se(value))

#data_barplot_melt_sum$ident <- factor(data_barplot_melt_sum$ident, levels = unique(data_barplot_melt_sum$ident))

#data_barplot_melt_sum$variable <- factor(data_barplot_melt_sum$variable, levels = unique(MarkerParkinsonsGenes))

#write.table(data_barplot_melt_sum, file = paste("Files/MajorCellTypes_ImmuneMarker_Genes_BarChartTable.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")







