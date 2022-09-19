#Workflow up to elbow plot to determine cutoff, and saving intermediate seurat
#object as .rds file in /n/scratch3/users/j/jap0606/batch1to8

#set.seed(100)

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

dir.create("Files/Resolution_Tests")

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls_5agebracket/FullIntegrationOSR_HealthyControls_5agebracket_MTG_Part1.rds")

SeuratObject <- FindNeighbors(SeuratObject, reduction = "harmony", dims = 1:50)

for(res in seq(from=0.5,to=1.4,by=0.1)){
   ResObject <- FindClusters(SeuratObject, resolution = res, algorithm = 4, method = "igraph")
   ResObject <- RunUMAP(ResObject, reduction = "harmony", dims = 1:50)

   cells_per_cluster_table <- as.data.frame(table(ResObject$seurat_clusters))

   colnames(cells_per_cluster_table) <- c("Cluster","Number of Cells")

   write.table(cells_per_cluster_table, file = paste("Files/Resolution_Tests/Cells_Per_Unassigned_Cluster_Table_Resolution",res,".tsv",sep=""), quote = FALSE, row.names = FALSE, sep = "\t")

   samples_per_cluster_table <- group_by(ResObject@meta.data, seurat_clusters) %>% summarise(Sample_Count = length(unique(sample_id)))

   colnames(samples_per_cluster_table) <- c("Cluster","Number of Samples")

   write.table(samples_per_cluster_table, file = paste("Files/Resolution_Tests/Samples_Per_Unassigned_Cluster_Table_Resolution",res,".tsv",sep=""), quote = FALSE, row.names = FALSE, sep = "\t")

   SeuratObject_UMAP_Clusters <- DimPlot(ResObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) + 
                                 theme(axis.text = element_text(size=8),
                                 axis.title = element_text(size = 12),
                                 legend.text = element_text(size = 8),
                                 title = element_text(size = 12),
                                 legend.key.size = unit(0.4,"cm"))

   ggsave(SeuratObject_UMAP_Clusters, filename = paste("Files/Resolution_Tests/SeuratObject_UMAP_Clusters_Resolution",res,".pdf",sep=""), device = "pdf", width = 6, height = 4, units = "in")
}

