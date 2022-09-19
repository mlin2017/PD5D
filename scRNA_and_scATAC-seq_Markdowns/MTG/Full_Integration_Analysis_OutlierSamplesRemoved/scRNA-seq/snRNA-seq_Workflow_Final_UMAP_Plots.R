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
library(Nebulosa)
library(scCustomize)

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegrationOSR_MTG_SeuratFinal.rds")

Idents(SeuratObject) <- "seurat_clusters"

SeuratObject <- subset(SeuratObject, subset = MajorCellTypes != "Unknown_Cluster_66")

SeuratObject <- RunUMAP(SeuratObject, reduction = "harmony", dims = 1:60)

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons_1", `2` = "Oligodendrocytes", `3` = "GLU_Neurons_2", `4` = "GLU_Neurons_3", `5` = "GLU_Neurons_4", `6` = "Astrocytes", `7` = "GLU_Neurons_5", `8` = "GLU_Neurons_6",`9` = "GLU_Neurons_7",`10` = "GABA_Neurons_1", `11` = "GLU_Neurons_8", `12` = "Microglia",`13` = "GABA_Neurons_2",`14` = "OPCs",`15` = "GABA_Neurons_3", `16`= "GABA_Neurons_4", `17`="GABA_Neurons_5", `18`="Endothelial_Cells_1", `19`="Endothelial_Cells_2", `21` = "GABA_Neurons_6", `22` = "GLU_Neurons_9", `23` = "GLU_Neurons_10", `24` = "GLU_Neurons_11", `25` = "GLU_Neurons_12", `26` = "GLU_Neurons_13",`27` = "GLU_Neurons_14", `30` = "GABA_Neurons_7", `32` = "GLU_Neurons_15", `33` = "GLU_Neurons_16",`34` = "GABA_Neurons_8", `35` = "GLU_Neurons_17", `36` = "GABA_Neurons_9", `38` = "GABA_Neurons_10", `41` = "GABA_Neurons_11", `42` = "GABA_Neurons_12",`43` = "GABA_Neurons_13", `44` = "GLU_Neurons_18", `45` = "GLU_Neurons_19", `47` = "TEMRA_T_Cells", `48` = "GABA_Neurons_14", `50` = "GABA_Neurons_15", `54` = "GLU_Neurons_20")

SeuratObject_UMAP_Clusters_nolegend <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) +
  theme(axis.text = element_text(size=8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 4),
              title = element_text(size = 12),
              legend.key.size = unit(0.4,"cm"),
              legend.position = "none")

ggsave(SeuratObject_UMAP_Clusters_nolegend, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_UMAP_Clusters_nolegendALT.pdf",sep=""), device = "pdf", width = 10.3, height = 10, units = "in")


SeuratObject_UMAP_Clusters_nolabel <- DimPlot(SeuratObject, reduction = "umap", label = FALSE, pt.size = 0.01, label.size=2.5, repel = TRUE) +
  theme(axis.text = element_text(size=7),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 6),
              title = element_text(size = 12))


ggsave(SeuratObject_UMAP_Clusters_nolabel, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_UMAP_Clusters_nolabelALT.pdf",sep=""), device = "pdf", width = 12.5, height = 10, units = "in")


SeuratObject_UMAP_Clusters_nolabel <- DimPlot(SeuratObject, reduction = "umap", label = FALSE, pt.size = 0.01, label.size=2.5, repel = TRUE) +
  theme(axis.text = element_text(size=7),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 5),
              title = element_text(size = 12))


ggsave(SeuratObject_UMAP_Clusters_nolabel, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_UMAP_Clusters_nolabelALT5.pdf",sep=""), device = "pdf", width = 12.5, height = 10, units = "in")


SeuratObject_UMAP_Clusters_nolegend_nolabel <- DimPlot(SeuratObject, reduction = "umap", label = FALSE, pt.size = 0.01, label.size=2.5, repel = TRUE) +
  theme(axis.text = element_text(size=7),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 6),
              title = element_text(size = 12),
              legend.key.size = unit(0.8,"cm"),
              legend.position = "none")

ggsave(SeuratObject_UMAP_Clusters_nolegend_nolabel, filename = paste("Figures/Subcluster_Assigned_",SeuratObject@project.name,"_UMAP_Clusters_nolegend_nolabelALT.pdf",sep=""), device = "pdf", width = 10.3, height = 10, units = "in")


