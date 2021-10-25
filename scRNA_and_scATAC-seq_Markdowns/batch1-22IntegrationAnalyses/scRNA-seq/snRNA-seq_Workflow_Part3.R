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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/batch1-22/Batch1-22_MTG_Part2.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons", `2` = "GLU_Neurons", `3` = "Oligodendrocytes", `4` = "GLU_Neurons", `5` = "Astrocytes", `6` = "GLU_Neurons", `7` = "GLU_Neurons", `8` = "GABA_Neurons",`9` = "GABA_Neurons",`10` = "GLU_Neurons", `11` = "Microglia",`12` = "GABA_Neurons",`13` = "Unknown_Cluster_13",`14` = "OPCs",`15` = "GLU_Neurons", `16`= "Cajal_Retzius_Cells", `17`="GLU_GABA_Neurons", `18`="GABA_Neurons", `19`="GABA_Neurons",  `20`="Endothelial", `21` = "Endothelial", `22` = "Unknown_Cluster_22", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GLU_Neurons",`27` = "GABA_Neurons",`28` = "Unknown_Cluster_28",`29` = "CD8+_T_Cells",`30` = "Unknown_Cluster_30")
                   

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




