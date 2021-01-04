# Load libraries
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(hdf5r)
set.seed(1234)

# read atacseq files : 1) count_mat, 2) metadata, 3) chrom_assay
# counts1 is a matrix with rows are genomic regions (ex. chr1:631141-631806 ) 
# and columns are samples (ex. AAACGAAAGAACCATA-2 ). the content in the matrix is 
# the nummber of TNFs in each cell in the specific genomic region.

counts <- Read10X_h5(filename = "./atac_merge/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = './atac_merge/fragments.tsv.gz'
)
metadata <- read.csv(
  file = "./atac_merge/singlecell.csv",
  header = TRUE,
  row.names = 1
)
# Create Seurat object for the atacseq data
atac0 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata)

# Preprocess atacseq data: 
preprocess_attac <- function(atac) {
  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
  
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg19"
  
  # add the gene information to the object
  Annotation(atac) <- annotations
  
  # compute nucleosome signal score per cell
  atac <- NucleosomeSignal(object = atac)
  
  # compute TSS enrichment score per cell
  atac <- TSSEnrichment(object = atac, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
  atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments
  
  atac$high.tss <- ifelse(atac$TSS.enrichment > 2, 'High', 'Low')
  TSSPlot(atac, group.by = 'high.tss') + NoLegend()
  
  
  atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  FragmentHistogram(object = atac, group.by = 'nucleosome_group')
  return(atac)
}
# Perform preprocessing steps outlined in the preprocess_attac function.
atac0<-preprocess_attac(atac0)
atac0 <- RunTFIDF(atac0)

# Find top binary features for a given assay based on total number of cells containing feature.
# Can specify a minimum cell count, or a lower percentile bound.
dim(atac0@assays$peaks@counts)
atac0 <- FindTopFeatures(atac0, min.cutoff = 'q0')
dim(atac0@assays$peaks@counts)

atac0 <- RunSVD(atac0, n = 50, reduction.name = "svd") # n =  50 by default, saved as lsi
atac0 <- RunLSI(atac0, n = 50, reduction.name = "lsi") # n =  50 by default, saved as lsi
atac0 <- RunUMAP(object = atac0, reduction = 'lsi', dims = 1:30)
atac0 <- FindNeighbors(object = atac0, reduction = 'lsi', dims = 1:30)
atac0 <- FindClusters(object = atac0, verbose = FALSE, algorithm = 3) # 3 is SLM algorithm
DimPlot(object = atac0,reduction = 'umap', label = TRUE) + NoLegend()

#summary above : 
# 1) find top 50 features,
# 2) do dimensionality reduction using svd (top 50) or lsi or lsi+umap
# 3) Find Neighbors and Clustering
# 4) Plots

atac0@reductions
saveRDS(atac0,"./atac_merge.rds")
###### 


atac0=readRDS("./atac_merge.rds")

# Compute counts per cell in gene body and promoter region. 
# add the gene activity matrix to the Seurat object as a new assay and normalize it
gene.activities <- GeneActivity(atac0)  
atac0[['RNA']] <- CreateAssayObject(counts = gene.activities)
saveRDS(atac0,"./atac_merge.rds")

#######
atac0=readRDS("./atac_merge.rds")

atac0$nCount_peaks
atac0$nCount_peaks
summary(atac0$nCount_peaks)
atac1 <- subset(atac0, subset = nCount_peaks >100)


DefaultAssay(atac1) <- 'RNA'
atac1 <- FindVariableFeatures(atac1)

atac1<- NormalizeData(
  object = atac1,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac1$nCount_RNA)
)


FeaturePlot(
  object = atac1,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  reduction = 'umap'
)

saveRDS(atac0,"./atac_merge.rds")
saveRDS(atac1,"./atac_merge_filtered.rds")
###############
atac1=readRDS("./atac_merge_filtered.rds")
rna=readRDS("single_cell_data/AllMB.rds")

rna@assays$RNA@counts
atac1@assays$RNA@counts
dim(rna@assays$RNA@counts)
dim(atac1@assays$RNA@counts)


# DefaultAssay(atac1) <- "peaks"
# VariableFeatures(atac1) <- names(which(Matrix::rowSums(atac1) > 100))
# atac1 <- RunLSI(atac1, n = 50, scale.max = NULL)
# atac1 <- RunUMAP(atac1, reduction = "lsi", dims = 1:50)
# p1 <- DimPlot(atac1, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
# p2 <- DimPlot(rna, label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
# p1 + p2


DefaultAssay(atac1)<"RNA"
transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = atac1,
  features = VariableFeatures(object = atac1), 
  reference.assay = "RNA",
  query.assay = "RNA",
  reduction = 'cca'
)

dim(transfer.anchors@anchors)
transfer.anchors@anchor.features
# rna$celltype <- Idents(rna)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = Idents(rna),
  weight.reduction = atac1[['lsi']],
  dims = 1:30
)

atac1 <- AddMetaData(object = atac1, metadata = predicted.labels)
saveRDS(atac1,"./atac_merge_filtered_matched_100.rds")
###############
rna=readRDS("single_cell_data/AllMB.rds")
atac1=readRDS("./atac_merge_filtered_matched_100.rds")


hist(atac1$prediction.score.max)
abline(v = 0.5, col = "red")

table(atac1$prediction.score.max > 0.1)
atac1.atac.filtered <- subset(atac1, subset = prediction.score.max > 0.1)
atac1.atac.filtered$predicted.id <- factor(atac1.atac.filtered$predicted.id, levels = levels(rna))  # to make 
p1 <- DimPlot(atac1.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(rna,group.by = "celltype",  label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()

p1 + p2




### Coembedding :

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
rna$celltype=Idents(rna)
genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac1[["lsi"]])

# this line adds the imputed data matrix to the atac object
atac1[["RNA"]] <- imputation
rna$tech="rna"
atac1$tech="atac"
coembed <- merge(x = rna, y = atac1)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)
coembed$nFeature_peaks

p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()


coembed$blacklist_region_fragments[is.na(coembed$blacklist_region_fragments)] <- 0
FeaturePlot(coembed, features = "blacklist_region_fragments", max.cutoff = 500)

##########
# change back to working with peaks instead of gene activities
DefaultAssay(atac1.atac.filtered) <- 'peaks'
Idents(atac1.atac.filtered)<-atac1.atac.filtered$predicted.id
da_peaks <- FindMarkers(
  object = atac1.atac.filtered,
  ident.1 = "dopamine_neuron_potentially_DA2",
  ident.2 = "Oligodendrocytes",
  min.pct = 0.01,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
saveRDS(da_peaks,"da_peaks.rds")
head(da_peaks)

plot1 <- VlnPlot(
  object = atac1.atac.filtered,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("dopamine_neuron_potentially_DA2","Oligodendrocytes")
)
plot2 <- FeaturePlot(
  object = atac1.atac.filtered,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

fc <- FoldChange(atac1.atac.filtered, ident.1 = "dopamine_neuron_potentially_DA2", ident.2 = "Oligodendrocytes")
head(fc)

open_DA2 <- rownames(da_peaks[da_peaks$avg_logFC > 0.5, ])
open_Oligodendrocytes <- rownames(da_peaks[da_peaks$avg_logFC < -0.5, ])

closest_genes_DA2 <- ClosestFeature(atac1.atac.filtered, regions = open_DA2)
closest_genes_Oligodendrocytes  <- ClosestFeature(atac1.atac.filtered, regions = open_Oligodendrocytes)

head(closest_genes_DA2)


head(closest_genes_Oligodendrocytes)

CoveragePlot(
  object = atac1.atac.filtered,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)

CoverageBrowser(
  object = atac1.atac.filtered,
  region = rownames(da_peaks)[1:10],
  extend.upstream = 40000,
  extend.downstream = 20000
)

library(shiny)
