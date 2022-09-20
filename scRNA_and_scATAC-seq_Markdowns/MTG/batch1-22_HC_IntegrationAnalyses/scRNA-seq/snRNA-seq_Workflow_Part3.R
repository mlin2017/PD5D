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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/batch1-22_HC/Batch1-22_MTG_HC_Part2.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons", `2` = "GLU_Neurons", `3` = "Oligodendrocytes", `4` = "GLU_Neurons", `5` = "GLU_Neurons", `6` = "GLU_Neurons", `7` = "GABA_Neurons", `8` = "Astrocytes",`9` = "GLU_Neurons",`10` = "GABA_Neurons", `11` = "GABA_Neurons",`12` = "GLU_Neurons",`13` = "Microglia",`14` = "Cajal_Retzius_Cells",`15` = "GLU_Neurons", `16`= "OPCs", `17`="Endothelial_Cells", `18`="Endothelial", `19`="GLU_Neurons",  `20`="GLU_Neurons", `21` = "GABA_Neurons", `22` = "GABA_Neurons", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "Oligo_GLU_Neurons", `26` = "GLU_Neurons",`27` = "GABA_Neurons",`28` = "GLU_Neurons")
                   
cells_per_cluster_table <- as.data.frame(table(SeuratObject@active.ident))

colnames(cells_per_cluster_table) <- c("Cell_Type","Number of Cells")

write.table(cells_per_cluster_table, file = "Files/Cells_Per_Cell_Type_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

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

"/n/scratch3/users/j/jap0606/batch1-22/Batch1to22_MTG_Part2.rds"


clusteridents <- rbind("/n/scratch3/users/j/jap0606/batch1-22/Batch1-22_MTG_Part2.rds",data.frame(as.vector(unique(SeuratObject@active.ident))))

write.table(clusteridents, "Files/ClusterIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

write.table(data.frame(as.vector(unique(SeuratObject@meta.data$sample_id))), "Files/SampleIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

SeuratObject_Case_Split_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", split.by = "case", label = TRUE, repel = TRUE, ncol = 1, label.size=2) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_Case_Split_UMAP_Clusters, filename = paste("Figures/Assigned_",SeuratObject@project.name,"_Case_Split_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 8, units = "in")

SeuratObject_Region_Split_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", split.by = "batch", label = TRUE, repel = TRUE, ncol = 1, label.size=2) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_Region_Split_UMAP_Clusters, filename = paste("Figures/Assigned_",SeuratObject@project.name,"_MTG_Region_Split_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 32, units = "in")

ggsave(SeuratObject_MTG_UMAP_Clusters, filename = paste("Figures/Assigned_",SeuratObject@project.name,"_MTG_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 4, units = "in")




