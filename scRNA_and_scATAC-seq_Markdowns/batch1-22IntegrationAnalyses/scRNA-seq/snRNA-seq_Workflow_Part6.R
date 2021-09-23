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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/batch1-22/Batch1to22_MTG_Part2.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons", `2` = "GLU_Neurons", `3` = "Oligodendrocytes", `4` = "GLU_Neurons", `5` = "Astrocytes", `6` = "GLU_Neurons", `7` = "GLU_Neurons", `8` = "GABA_Neurons",`9` = "GABA_Neurons",`10` = "GLU_Neurons", `11` = "Microglia",`12` = "GABA_Neurons",`13` = "Unknown_Cluster_13",`14` = "OPCs",`15` = "GLU_Neurons", `16`= "Cajal_Retzius_Cells", `17`="GLU_GABA_Neurons", `18`="GABA_Neurons", `19`="GABA_Neurons",  `20`="Endothelial", `21` = "Endothelial", `22` = "Unknown_Cluster_22", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GLU_Neurons",`27` = "GABA_Neurons",`28` = "Unknown_Cluster_28",`29` = "CD8+_T_Cells",`30` = "Unknown_Cluster_30")

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

write.table(data_barplot_melt_sum, file = paste("Files/Parkisons_Mendelian_Genes_BarChartTable.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
