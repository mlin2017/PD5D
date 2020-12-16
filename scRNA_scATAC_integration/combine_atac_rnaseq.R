library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)

'atac_merge'
#####
counts1<- Read10X_h5(filename = "./atac_merge/filtered_peak_bc_matrix.h5")
#####
metadata1 <- read.csv(
  file = "./atac_merge/singlecell.csv",
  header = TRUE,
  row.names = 1
)

#####
chrom_assay1 <- CreateChromatinAssay(
  counts = counts1,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = './atac_merge/fragments.tsv.gz'
)
#####
#####
pbmc1 <- CreateSeuratObject(
  counts = chrom_assay1,
  assay = "peaks",
  meta.data = metadata1)
#####
preprocess_attac <- function(pbmc) {
  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
  
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg19"
  
  # add the gene information to the object
  Annotation(pbmc) <- annotations
  
  # compute nucleosome signal score per cell
  pbmc <- NucleosomeSignal(object = pbmc)
  
  # compute TSS enrichment score per cell
  pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
  pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
  
  pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
  TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
  
  
  pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
  return(pbmc)
}
pbmc1<-preprocess_attac(pbmc1)
pbmc1 <- RunTFIDF(pbmc1)
pbmc1 <- FindTopFeatures(pbmc1, min.cutoff = 'q0')
pbmc1 <- RunSVD(pbmc1)
pbmc1 <- RunUMAP(object = pbmc1, reduction = 'lsi', dims = 2:30)
pbmc1 <- FindNeighbors(object = pbmc1, reduction = 'lsi', dims = 2:30)
pbmc1 <- FindClusters(object = pbmc1, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc1, label = TRUE) + NoLegend()

saveRDS(pbmc1,"./pbmc_atac_merge.rds")


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)

pbmc=readRDS("./pbmc_atac_merge.rds")

gene.activities <- GeneActivity(pbmc)

saveRDS(pbmc,"./pbmc_atac_merge_orig.rds")
saveRDS(gene.activities,"./gene.activities_pbmc_atac_merge.rds")


# pbmc<-pbmc1
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
saveRDS(pbmc,"./pbmc_atac_merge.rds")

summary(pbmc$nCount_peaks)
pbmc <- subset(pbmc, subset = nCount_peaks > 300)


DefaultAssay(pbmc) <- 'RNA'
pbmc <- FindVariableFeatures(pbmc)

pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)


FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)


saveRDS(pbmc,"./pbmc_atac_merge.rds")


pbmc=readRDS("./pbmc_atac_merge.rds")
AllMB=readRDS("single_cell_data/AllMB.rds")

AllMB@assays$RNA@counts
dim(pbmc@assays$RNA@counts)

transfer.anchors <- FindTransferAnchors(
  reference = AllMB,
  query = pbmc,
  features = VariableFeatures(object = pbmc), 
  reference.assay = "RNA",
  query.assay = "RNA",
  reduction = 'cca'
)

# AllMB$celltype <- Idents(AllMB)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = Idents(AllMB),
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
saveRDS(pbmc,"./pbmc_atac_merge_matched.rds")
###############
AllMB=readRDS("single_cell_data/AllMB.rds")
pbmc=readRDS("./pbmc_atac_merge_matched.rds")


table(pbmc$prediction.score.max > 0.5)
pbmc.atac.filtered <- subset(pbmc, subset = prediction.score.max > 0.5)
pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, levels = levels(AllMB))  # to make 
p1 <- DimPlot(pbmc.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(AllMB, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
p1 + p2




### Coembedding :

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(AllMB)
refdata <- GetAssayData(AllMB, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc[["lsi"]])

# this line adds the imputed data matrix to the pbmc object
pbmc[["RNA"]] <- imputation
AllMB$tech="rna"
pbmc$tech="atac"
coembed <- merge(x = AllMB, y = pbmc)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()


# AllMB$celltype <- Idents(AllMB)
# plot1 <- DimPlot(
#   object = AllMB,
#   group.by = 'celltype',
#   label = TRUE,
#   repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
# 
# plot2 <- DimPlot(
#   object = pbmc,
#   group.by = 'predicted.id',
#   label = TRUE,
#   repel = FALSE) + ggtitle('scATAC-seq')
# 
# plot1 + plot2
# 
# pbmc$predicted.id
