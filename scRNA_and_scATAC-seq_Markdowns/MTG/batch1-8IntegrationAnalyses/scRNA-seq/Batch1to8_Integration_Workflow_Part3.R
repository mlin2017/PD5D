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
library(stringr)
library(reshape2)
library(sciplot)

Batch1to8_MTG <- readRDS("/n/scratch3/users/j/jap0606/batch1to8/Batch1to8_MTG_Part2.rds")

Batch1to8_MTG <- RenameIdents(Batch1to8_MTG, `1` = "GLU_Neurons_1", `2` = "Oligodendrocytes", `3` = "GLU_Neurons_2", `4` = "GLU_Neurons_3", `5` = "Astrocytes", `6` = "GABA_Neurons_1", `7` = "GLU_Neurons_4", `8` = "GLU_Neurons_5",`9` = "GABA_Neurons_2",`10` = "Microglia", `11` = "GLU_Neurons_6",`12` = "GLU_Neurons_7",`13` = "GABA_Neurons_3",`14` = "OPCs",`15` = "GLU_Neurons_8", `16`= "GLU_Neurons_9", `17`="GABA_Neurons_4", `18`="GABA_Neurons_5", `19`="Endothelial_Cells_1",  `20`="Endothelial_Cells_2", `21` = "GLU_Neurons_10", `22` = "GLU_Neurons_11", `23` = "GLU_Neurons_12", `24` = "GLU_Neurons_13", `25` = "GLU_Neurons_14", `26` = "GABA_Neurons_6",`27` = "Cajal_Retzius_Cells",`28` = "GLU_Neurons_15",`29` = "Unknown_Cluster_29",`30` = "Unknown_Cluster_30")
                   

Batch1to8_MTG_UMAP_Clusters <- DimPlot(Batch1to8_MTG, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

"/n/scratch3/users/j/jap0606/batch1to8/Batch1to8_MTG_Part2.rds"


clusteridents <- rbind("/n/scratch3/users/j/jap0606/batch1to8/Batch1to8_MTG_Part2.rds",data.frame(as.vector(unique(Batch1to8_MTG@active.ident))))

write.table(clusteridents, "Files/ClusterIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

write.table(data.frame(as.vector(unique(Batch1to8_MTG@meta.data$sample_id))), "Files/SampleIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

Batch1to8_MTG_Case_Split_UMAP_Clusters <- DimPlot(Batch1to8_MTG, reduction = "umap", split.by = "case", label = TRUE, repel = TRUE, ncol = 1, label.size=2) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(Batch1to8_MTG_Case_Split_UMAP_Clusters, filename = paste("Figures/Assigned_",Batch1to8_MTG@project.name,"_Case_Split_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 8, units = "in")

Batch1to8_MTG_Region_Split_UMAP_Clusters <- DimPlot(Batch1to8_MTG, reduction = "umap", split.by = "batch", label = TRUE, repel = TRUE, ncol = 1, label.size=2) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(Batch1to8_MTG_Region_Split_UMAP_Clusters, filename = paste("Figures/Assigned_",Batch1to8_MTG@project.name,"_MTG_Region_Split_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 32, units = "in")

ggsave(Batch1to8_MTG_UMAP_Clusters, filename = paste("Figures/Assigned_",Batch1to8_MTG@project.name,"_MTG_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 4, units = "in")

saveRDS(Batch1to8_MTG,"/n/scratch3/users/j/jap0606/batch1to8/Batch1to8_MTG_ClustersAssigned.rds")


