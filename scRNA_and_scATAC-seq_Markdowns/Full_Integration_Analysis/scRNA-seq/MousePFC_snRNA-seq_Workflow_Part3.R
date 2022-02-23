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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/Mouse_PFC_Seurat/Mouse_PFC_Part3.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "GABA_Neurons", `2` = "Astrocytes", `3` = "GABA_Neurons", `4` = "GLU_Neurons", `5` = "GLU_Neurons", `6` = "GLU_Neurons", `7` = "GLU_Neurons", `8` = "GLU_Neurons",`9` = "Oligodendrocytes",`11` = "GLU_Neurons",`12` = "GABA_Neurons",`13` = "GABA_Neurons",`14` = "OPCs",`15` = "GABA_Neurons", `16`= "GABA_Neurons", `17`="GLU_Cajal_Retzius", `18`="Endothelial_Cells", `19`="GLU_Neurons",  `20`="GLU_Neurons", `21` = "GABA_Neurons", `22` = "Oligo_Endo", `23` = "Microglia", `24` = "Ependymal_Cells", `26` = "GLU_Neurons",`27` = "GLU_Neurons")


CellType_SplitCase_ProportionTable <- as.data.frame(prop.table(table(Idents(SeuratObject), SeuratObject$case), margin = 2))

write.table(CellType_SplitCase_ProportionTable, file = "Files/CellType_SplitCase_ProportionTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellType_SplitSample_ProportionTable <- as.data.frame(prop.table(table(Idents(SeuratObject), SeuratObject$sample_ID), margin = 2))

write.table(CellType_SplitSample_ProportionTable, file = "Files/CellType_SplitSample_ProportionTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellType_SplitCase_NumbersTable <- as.data.frame(table(Idents(SeuratObject), SeuratObject$case))

write.table(CellType_SplitCase_NumbersTable, file = "Files/CellType_SplitCase_NumbersTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellType_SplitSample_NumbersTable <- as.data.frame(table(Idents(SeuratObject), SeuratObject$sample_ID))

write.table(CellType_SplitSample_NumbersTable, file = "Files/CellType_SplitSample_NumbersTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

                   
SeuratObject_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

clusteridents <- rbind("/n/scratch3/users/j/jap0606/Mouse_PFC_Seurat/Mouse_PFC_Part3.rds",data.frame(as.vector(unique(SeuratObject@active.ident))))

write.table(clusteridents, "Files/ClusterIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

write.table(data.frame(as.vector(unique(SeuratObject@meta.data$sample_ID))), "Files/SampleIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

SeuratObject_Case_Split_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", split.by = "case", label = TRUE, repel = TRUE, ncol = 1, label.size=2) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_Case_Split_UMAP_Clusters, filename = paste("Figures/Assigned_",SeuratObject@project.name,"_Case_Split_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 8, units = "in")

ggsave(SeuratObject_UMAP_Clusters, filename = paste("Figures/Assigned_",SeuratObject@project.name,"_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 4, units = "in")

#saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/Mouse_PFC_Seurat/Mouse_PFC_Part3.rds")


