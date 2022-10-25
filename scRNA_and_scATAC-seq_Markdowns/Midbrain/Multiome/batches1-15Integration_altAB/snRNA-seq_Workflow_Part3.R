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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_PostClusterQCUMAP.rds")

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons", `2` = "Oligodendrocytes", `3` = "GLU_Neurons", `4` = "GLU_Neurons", `5` = "GLU_Neurons", `6` = "Astrocytes", `7` = "GLU_Neurons", `8` = "GLU_Neurons",`9` = "GLU_Neurons",`10` = "GABA_Neurons", `11` = "GLU_Neurons", `12` = "Microglia",`13` = "GABA_Neurons",`14` = "OPCs",`15` = "GABA_Neurons", `16`= "GABA_Neurons", `17`="GABA_Neurons", `18`="Endothelial_Cells", `19`="Endothelial_Cells", `21` = "GABA_Neurons", `22` = "GLU_Neurons", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GLU_Neurons",`27` = "GLU_Neurons", `30` = "GABA_Neurons", `32` = "GLU_Neurons", `33` = "GLU_Neurons",`34` = "GABA_Neurons", `35` = "GLU_Neurons", `36` = "GABA_Neurons", `38` = "GABA_Neurons", `41` = "GABA_Neurons", `42` = "GABA_Neurons",`43` = "GABA_Neurons", `44` = "GLU_Neurons", `45` = "GLU_Neurons", `47` = "TEMRA_T_Cells", `48` = "GABA_Neurons", `50` = "GABA_Neurons", `54` = "GLU_Neurons", `66` = "Unknown_Cluster_66")


CellType_SplitCase_ProportionTable <- as.data.frame(prop.table(table(Idents(SeuratObject), SeuratObject$case), margin = 2))

write.table(CellType_SplitCase_ProportionTable, file = "Files/CellType_SplitCase_ProportionTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellType_SplitSample_ProportionTable <- as.data.frame(prop.table(table(Idents(SeuratObject), SeuratObject$sample_id), margin = 2))

write.table(CellType_SplitSample_ProportionTable, file = "Files/CellType_SplitSample_ProportionTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellType_SplitCase_NumbersTable <- as.data.frame(table(Idents(SeuratObject), SeuratObject$case))

write.table(CellType_SplitCase_NumbersTable, file = "Files/CellType_SplitCase_NumbersTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellType_SplitSample_NumbersTable <- as.data.frame(table(Idents(SeuratObject), SeuratObject$sample_id))

write.table(CellType_SplitSample_NumbersTable, file = "Files/CellType_SplitSample_NumbersTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

                   
SeuratObject_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))


clusteridents <- rbind("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_PostClusterQCUMAP.rds",data.frame(as.vector(unique(SeuratObject@active.ident))))


SeuratObject_UMAP_Clusters_Individual <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) +
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))


write.table(clusteridents, "Files/ClusterIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

write.table(data.frame(as.vector(unique(SeuratObject@meta.data$sample_id))), "Files/SampleIdents.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

SeuratObject_Case_Split_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", split.by = "case", label = TRUE, repel = TRUE, ncol = 1, label.size=2) + 
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 8),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_Case_Split_UMAP_Clusters, filename = paste("Figures/Assigned_",SeuratObject@project.name,"_Case_Split_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 8, units = "in")

ggsave(SeuratObject_UMAP_Clusters, filename = paste("Figures/Assigned_",SeuratObject@project.name,"_UMAP_Clusters.pdf",sep=""), device = "pdf", width = 6, height = 4, units = "in")

celltable <- data.frame(table(Idents(SeuratObject), SeuratObject$case))

dcast(data = celltable,formula = Var1~Var2,fun.aggregate = sum,value.var = "Freq")

write.table(clusteridents, "Files/CellType_CaseSplit_FrequencyTable.txt", col.names =FALSE, quote = FALSE, row.names = FALSE)

cluster_assignment_table <- as.data.frame(cbind(rownames(SeuratObject@meta.data),as.vector(SeuratObject@active.ident)))

write.table(cluster_assignment_table, file = "Files/CellBarcodes_to_MajorCellType_Assignment_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

