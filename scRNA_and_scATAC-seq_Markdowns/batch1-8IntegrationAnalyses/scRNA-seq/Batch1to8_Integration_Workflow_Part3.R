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

Batch1to8_MTG <- readRDS("/n/scratch3/users/j/jap0606/batch1to8/Batch1to8_MTG_Part2.rds")

Batch1to8_MTG <- RenameIdents(Batch1to8_MTG, `0` = "GLU_Neurons", `1` = "GLU_Neurons", `2` = "Oligodendrocytes",
                      `3` = "GLU_Neurons", `4` = "Oligodendrocytes", `5` = "GABA_Neurons",
                      `6` = "GLU_Neurons", `7` = "Astrocytes", `8` = "GLU_Neurons",`9` = "GLU_Neurons",
                      `10` = "GABA_Neurons", `11` = "Microglia",`12` = "GABA_Neurons",
                      `13` = "GLU_Neurons",`14` = "GABA_Neurons",
                      `15` = "OPCs", `16`= "GLU_Neurons", `17`="Endothelial", `18`="Endothelial", `19`="GLU_Neurons",  `20`="GABA_Neurons", `21` = "Oligo_GLU_Neurons", `22` = "GLU_Neurons", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GLU_Neurons",`27` = "GABA_Neurons",`28` = "GLU_Neurons",`29` = "GLU_Neurons",`30` = "Oligo_GABA_Neurons",`31` = "GABA_Neurons",`32` = "GLU_Neurons",`33` = "GABA_Neurons",`34` = "Astrocytes", `35` = "GLU_Neurons")

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




