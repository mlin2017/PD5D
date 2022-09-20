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
library(tidyr)
library(scDblFinder)

dir.create("Figures/Cluster_QC")

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration/FullIntegration_MTG_Part2.rds")

SeuratObject.sce <- as.SingleCellExperiment(SeuratObject)
SeuratObject.sce <- scDblFinder(SeuratObject.sce, samples = "sample_id",clusters = "seurat_clusters")

barcode_doublet_Table <- as.data.frame(cbind(as.vector(rownames(SeuratObject.sce@colData)),as.vector(SeuratObject.sce@colData$scDblFinder.class)))

SeuratObject@meta.data$class <- barcode_doublet_Table$V2

Doublet_UMAP <- DimPlot(SeuratObject, label = FALSE, repel = TRUE, pt.size = 0, label.size = 2.5, group.by = "class") + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(Doublet_UMAP, filename = "Figures/Cluster_QC/Unassigned_Doublet_UMAPclusters_scRNA_seq.pdf", device = "pdf", width = 6, height = 4, units = "in")

singlet_barcodes <- rownames(SeuratObject@meta.data)[as.vector(SeuratObject@meta.data$class) %in% "singlet"]

Doublets_Removed_UMAP <- DimPlot(SeuratObject, label = TRUE, repel = TRUE, pt.size = 0, label.size = 2.5, cells = singlet_barcodes) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(Doublets_Removed_UMAP, filename = "Figures/Cluster_QC/Unassigned_DoubletsRemoved_UMAPclusters_scRNA_seq.pdf", device = "pdf", width = 6, height = 4, units = "in")

metadata <- SeuratObject@meta.data

doublet_counts <- metadata %>% count(seurat_clusters, class)

cor_doublet_counts <- data.frame(matrix(ncol = 3, nrow = 0))

colnames(cor_doublet_counts) <- colnames(doublet_counts)

for (i in unique(doublet_counts$seurat_clusters)) {
  temptable <- doublet_counts[doublet_counts$seurat_clusters %in% i,]
  if (length(temptable$seurat_clusters) == 2){
    cor_doublet_counts <- rbind(cor_doublet_counts, temptable)
  }
  else {
    cell_class <- temptable$class[1]
    if (cell_class == "singlet") {
      temptable2 <- temptable
      temptable2[2] <- "doublet"
      temptable2[3] <- 0
      cor_doublet_counts <- rbind(cor_doublet_counts, temptable, temptable2)
    }
    else if (cell_class == "doublet") {
      temptable2 <- temptable
      temptable2[2] <- "doublet"
      temptable2[3] <- 0
      cor_doublet_counts <- rbind(cor_doublet_counts, temptable, temptable2)
  }
 }
}

doublet_counts <- cor_doublet_counts

doublet_percentage <- doublet_counts %>% group_by(seurat_clusters) %>% mutate(doublet_percent = n[which(class == "doublet")]/sum(n)*100, total_cells = sum(n))

write.table(doublet_percentage, file = "Files/Doublet_Percentage_Per_Cluster_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

doublet_percentage <- doublet_percentage[doublet_percentage$class == "doublet",]

cell_threshold_clusters <- as.numeric(doublet_percentage$seurat_clusters[doublet_percentage$total_cells >= 100])

SeuratObjectSinglets <- subset(SeuratObject, subset = class == "singlet" & seurat_clusters %in% cell_threshold_clusters)

singlet_barcodes_df <- as.data.frame(rownames(SeuratObjectSinglets@meta.data))

write.table(singlet_barcodes_df, file = "Files/Singlet_Barcodes_Dataframe.tsv", quote = FALSE, row.names = FALSE, sep = "\t")



