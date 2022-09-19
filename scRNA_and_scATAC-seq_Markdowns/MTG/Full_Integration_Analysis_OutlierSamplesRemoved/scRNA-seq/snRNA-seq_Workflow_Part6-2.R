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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration/FullIntegrationOSR_MTG_PostClusterQCUMAP.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons_1", `2` = "Oligodendrocytes", `3` = "GLU_Neurons_2", `4` = "GLU_Neurons_3", `5` = "GLU_Neurons_4", `6` = "Astrocytes", `7` = "GLU_Neurons_5", `8` = "GLU_Neurons_6",`9` = "GLU_Neurons_7",`10` = "GABA_Neurons_1", `11` = "GLU_Neurons_8", `12` = "Microglia",`13` = "GABA_Neurons_2",`14` = "OPCs",`15` = "GABA_Neurons_3", `16`= "GABA_Neurons_4", `17`="GABA_Neurons_5", `18`="Endothelial_Cells_1", `19`="Endothelial_Cells_2", `21` = "GABA_Neurons_6", `22` = "GLU_Neurons_9", `23` = "GLU_Neurons_10", `24` = "GLU_Neurons_11", `25` = "GLU_Neurons_12", `26` = "GLU_Neurons_13",`27` = "GLU_Neurons_14", `30` = "GABA_Neurons_7", `32` = "GLU_Neurons_15", `33` = "GLU_Neurons_16",`34` = "GABA_Neurons_8", `35` = "GLU_Neurons_17", `36` = "GABA_Neurons_9", `38` = "GABA_Neurons_10", `41` = "GABA_Neurons_11", `42` = "GABA_Neurons_12",`43` = "GABA_Neurons_13", `44` = "GLU_Neurons_18", `45` = "GLU_Neurons_19", `47` = "TEMRA_T_Cells", `48` = "GABA_Neurons_14", `50` = "GABA_Neurons_15", `54` = "GLU_Neurons_20", `66` = "Unknown_Cluster_66")

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

write.table(data_barplot_melt_sum, file = paste("Files/Clusters_Parkisons_Mendelian_Genes_BarChartTable.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

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

write.table(data_barplot_melt_sum, file = paste("Files/Clusters_Monika_Genes_BarChartTable.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

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

write.table(data_barplot_melt_sum, file = paste("Files/Clusters_Bea_Genes_BarChartTable.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

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







