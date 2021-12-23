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

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/batch1-22/Batch1-22_MTG_Part2.rds")

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

ggsave(Doublet_UMAP, filename = "Figures/Unassigned_Doublet_UMAPclusters_scRNA_seq.pdf", device = "pdf", width = 6, height = 4, units = "in")

singlet_barcodes <- rownames(SeuratObject@meta.data)[as.vector(SeuratObject@meta.data$class) %in% "singlet"]

Doublets_Removed_UMAP <- DimPlot(SeuratObject, label = TRUE, repel = TRUE, pt.size = 0, label.size = 2.5, cells = singlet_barcodes) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(Doublets_Removed_UMAP, filename = "Figures/Unassigned_DoubletsRemoved_UMAPclusters_scRNA_seq.pdf", device = "pdf", width = 6, height = 4, units = "in")

metadata <- SeuratObject@meta.data

doublet_counts <- metadata %>% count(seurat_clusters, class)

doublet_percentage <- doublet_counts %>% group_by(seurat_clusters) %>% mutate(doublet_percent = n[which(class == "doublet")]/sum(n)*100, total_cells = sum(n))

write.table(doublet_percentage, file = "Files/Doublet_Percentage_Per_Cluster_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

doublet_percentage <- doublet_percentage[doublet_percentage$class == "doublet",]
VlnPlot(SeuratObject, features = c("nCount_RNA"), ncol = 3, split.by = "seurat_clusters")

ncount_per_cluster <- VlnPlot(SeuratObject, features = c("nCount_RNA"), split.by = "seurat_clusters")

ggsave(ncount_per_cluster, filename = "Figures/nCount_per_Cluster_ViolinPlot.pdf", device = "pdf", width = 8, height = 6, units = "in")
mean_ncounts <- metadata %>% group_by(seurat_clusters) %>% summarise(mean = mean(nCount_RNA), SE = se(nCount_RNA))

median_ncounts <- metadata %>% group_by(seurat_clusters) %>% summarise(median = median(nCount_RNA), SE = se(nCount_RNA))
VlnPlot(SeuratObject, features = c("nFeature_RNA"), ncol = 3, split.by = "seurat_clusters")

nfeature_per_cluster <- VlnPlot(SeuratObject, features = c("nFeature_RNA"), split.by = "seurat_clusters")

ggsave(nfeature_per_cluster, filename = "Figures/nFeature_per_Cluster_ViolinPlot.pdf", device = "pdf", width = 8, height = 6, units = "in")



