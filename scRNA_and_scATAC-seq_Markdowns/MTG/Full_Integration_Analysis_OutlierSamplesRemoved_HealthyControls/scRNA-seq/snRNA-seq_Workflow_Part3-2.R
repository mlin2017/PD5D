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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls/FullIntegrationOSR_HealthyControls_MTG_PostClusterQCUMAP.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons_1", `2` = "Oligodendrocytes", `3` = "GLU_Neurons_2", `4` = "GLU_Neurons_3", `5` = "GLU_Neurons_4", `6` = "Astrocytes", `7` = "GLU_Neurons_5", `8` = "GLU_Neurons_6",`9` = "GLU_Neurons_7",`10` = "GABA_Neurons_1", `11` = "Microglia", `12` = "GLU_Neurons_8",`13` = "GLU_Neurons_9",`14` = "GABA_Neurons_2",`15` = "Endothelial_Cells_1", `16`= "OPCs", `17`="GABA_Neurons_3", `18` = "GLU_Neurons_10", `19`="GABA_Neurons_4", `20`="Endothelial_Cells_2", `21` = "GABA_Neurons_5", `23` = "GLU_Neurons_11", `24` = "GABA_Neurons_6", `26` = "GABA_Neurons_7",`27` = "GLU_Neurons_12", `28` = "GLU_Neurons_13", `29` = "GLU_Neurons_14", `31` = "GLU_Neurons_15", `32` = "GABA_Neurons_8",`34` = "GABA_Neurons_9", `35` = "GABA_Neurons_10", `37` = "GABA_Neurons_11")


CellType_SplitCase_ProportionTable <- as.data.frame(prop.table(table(Idents(SeuratObject))))

write.table(CellType_SplitCase_ProportionTable, file = "Files/Subcluster_ProportionTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellType_SplitSample_ProportionTable <- as.data.frame(prop.table(table(Idents(SeuratObject), SeuratObject$sample_id), margin = 2))

write.table(CellType_SplitSample_ProportionTable, file = "Files/Subcluster_SplitSample_ProportionTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellType_SplitCase_NumbersTable <- as.data.frame(table(Idents(SeuratObject)))

write.table(CellType_SplitCase_NumbersTable, file = "Files/Subcluster_NumbersTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellType_SplitSample_NumbersTable <- as.data.frame(table(Idents(SeuratObject), SeuratObject$sample_id))

write.table(CellType_SplitSample_NumbersTable, file = "Files/Subcluster_SplitSample_NumbersTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

                   
SeuratObject_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))


SeuratObject_UMAP_Clusters_nolabel <- DimPlot(SeuratObject, reduction = "umap", label = FALSE, pt.size = 0.01, label.size=2, repel = TRUE) +
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 4),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_UMAP_Clusters_nolabel, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_UMAP_Clusters_nolabels.pdf",sep=""), device = "pdf", width = 12, height = 4, units = "in")

SeuratObject_UMAP_Clusters_nolegend <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2, repel = TRUE) +
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 4),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"),
              legend.position = "none")

ggsave(SeuratObject_UMAP_Clusters_nolegend, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_UMAP_Clusters_nolegend.pdf",sep=""), device = "pdf", width = 10, height = 10, units = "in")


SeuratObject_UMAP_Clusters_BIG <- DimPlot(SeuratObject, reduction = "umap", label = FALSE, pt.size = 0.01, label.size=2, repel = TRUE) +
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 4),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))


ggsave(SeuratObject_UMAP_Clusters_BIG, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_UMAP_Clusters_BIG.pdf",sep=""), device = "pdf", width = 16, height = 10, units = "in")


clusteridents <- rbind("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls_5agebracket/FullIntegrationOSR_HealthyControls_5agebracket_MTG_PostClusterQCUMAP.rds",data.frame(as.vector(unique(SeuratObject@active.ident))))

write.table(clusteridents, "Files/SubclusterIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

write.table(data.frame(as.vector(unique(SeuratObject@meta.data$sample_id))), "Files/SampleIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

#SeuratObject_Case_Split_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", split.by = "case", label = TRUE, repel = TRUE, ncol = 1, label.size=2) + 
#  theme(axis.text = element_text(size=8),
#              axis.title = element_text(size = 12),
#              legend.text = element_text(size = 4),
#              title = element_text(size = 12),
#              legend.key.size = unit(0.4,"cm"))

#ggsave(SeuratObject_Case_Split_UMAP_Clusters, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_Case_Split_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 8, units = "in")

ggsave(SeuratObject_UMAP_Clusters, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 4, units = "in")

#celltable <- data.frame(table(Idents(SeuratObject), SeuratObject$case))

#dcast(data = celltable,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")

#write.table(clusteridents, "Files/Cluster_CaseSplit_FrequencyTable.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

cluster_assignment_table <- as.data.frame(cbind(rownames(SeuratObject@meta.data),as.vector(SeuratObject@active.ident)))

write.table(cluster_assignment_table, file = "Files/CellBarcodes_to_CellCluster_Assignment_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
#saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/Mouse_PFC_Seurat/Mouse_PFC_Part3.rds")


