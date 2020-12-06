install.packages("hdf5r")
install.packages("Signac")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools", "S4Vectors", "Biobase", "BiocGenerics", "Biostrings", "ggbio", "biovizBase", "AnnotationFilter"))
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)

'atac_merge'
#####
counts1 <- Read10X_h5(filename = "./attac_seq_data/B0085/filtered_peak_bc_matrix.h5")
counts2 <- Read10X_h5(filename = "./attac_seq_data/B0085/filtered_peak_bc_matrix.h5")
#####
metadata1 <- read.csv(
  file = "./attac_seq_data/B0085/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata2 <- read.csv(
  file = "./attac_seq_data/B0085/singlecell.csv",
  header = TRUE,
  row.names = 1
)
#####
chrom_assay1 <- CreateChromatinAssay(
  counts = counts1,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = './attac_seq_data/B0085/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
chrom_assay2 <- CreateChromatinAssay(
  counts = counts2,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = './attac_seq_data/B0085/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

# 
# counts <- Read10X_h5(filename = "./atac_merge/filtered_peak_bc_matrix.h5")
# 
# metadata <- read.csv(
#   file = "./atac_merge/singlecell.csv",
#   header = TRUE,
#   row.names = 1
# )
# 
# counts
# 
# chrom_assay <- CreateChromatinAssay(
#   counts = counts,
#   sep = c(":", "-"),
#   genome = 'hg19',
#   fragments = './atac_merge/fragments.tsv.gz',
#   min.cells = 1,  #10
#   min.features = 1  #200
# )
# chrom_assay


pbmc1 <- CreateSeuratObject(
  counts = chrom_assay1,
  assay = "peaks",
  meta.data = metadata1)

pbmc2 <- CreateSeuratObject(
  counts = chrom_assay2,
  assay = "peaks",
  meta.data = metadata2)

pbmc1$dataset="A"
pbmc2$dataset="B"

combined <- merge(
  x = pbmc1,
  y = list(pbmc2),
  add.cell.ids = c("A", "B")
)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)


pbmc.peaks

pbmc[['peaks']]


granges(pbmc)

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

VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


pbmc_orig<-pbmc
# pbmc <- subset(
#   x = pbmc,
#   subset = peak_region_fragments > 3000 &
#     peak_region_fragments < 20000 &
#     pct_reads_in_peaks > 15 &
#     blacklist_ratio < 0.05 &
#     nucleosome_signal < 4 &
#     TSS.enrichment > 2
# )
pbmc

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

DepthCor(pbmc)


pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()
gene.activities <- GeneActivity(pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)
DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('SNCA', 'GATA2', 'FOXA2', 'TH', 'TREM1', 'GFAP'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

head(pbmc$orig.ident)
pbmc$nCount_peaks
pbmc$nFeature_peaks
pbmc$TSS_fragments
pbmc$DNase_sensitive_region_fragments
pbmc$enhancer_region_fragments
pbmc$promoter_region_fragments
pbmc$on_target_fragments
pbmc$blacklist_region_fragments
pbmc$peak_region_fragments
pbmc$peak_region_cutsites
pbmc$passed_filters
pbmc$duplicate
pbmc$cell_id
pbmc$is__cell_barcode
pbmc$nucleosome_signal
pbmc$nucleosome_percentile
pbmc$TSS.enrichment
pbmc$TSS.percentile
pbmc$pct_reads_in_peaks
pbmc$blacklist_ratio
pbmc$high.tss
pbmc$nucleosome_group
pbmc$peaks_snn_res.0.8
pbmc$seurat_clusters
pbmc$nCount_RNA
pbmc$nFeature_RNA
assay(pbmc)
pbmc.markers


saveRDS(pbmc,'pbmc_B0085.rds')
saveRDS(pbmc,'pbmc_merge.rds')
save(pbmc,file="pbmc_merge.Robj")
pbmc<-readRDS('pbmc_merge.rds')
# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "1",
  ident.2 = "2",
  min.pct = 0.02,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)

plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("1","2")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2
########
plot1 <- VlnPlot(
  object = pbmc,
  features = 'SNCA',
  pt.size = 0.1,
  idents = c("6","7")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = 'SNCA',
  pt.size = 0.1
)

plot1 | plot2

########

fc <- FoldChange(pbmc, ident.1 = "1", ident.2 = "2")
head(fc)

head(da_peaks)
open_cd4naive <- rownames(da_peaks[da_peaks$avg_logFC > 0.5, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_logFC < -0.5, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)
head(closest_genes_cd4naive)

########################################################################
########################################################################
########################################################################
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("/home/stuartt/github/chrom/vignette_data/pbmc_10k_v3.rds")

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2


pbmc <- subset(pbmc, idents = 14, invert = TRUE)
pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'CD14 Mono',
  '1' = 'CD4 Memory',
  '2' = 'CD8 Effector',
  '3' = 'CD4 Naive',
  '4' = 'CD14 Mono',
  '5' = 'DN T',
  '6' = 'CD8 Naive',
  '7' = 'NK CD56Dim',
  '8' = 'pre-B',
  '9' = 'CD16 Mono',
  '10' = 'pro-B',
  '11' = 'DC',
  '12' = 'NK CD56bright',
  '13' = 'pDC'
)

# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14 Mono",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)


plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Memory","CD14 Mono")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

fc <- FoldChange(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14 Mono")
head(fc)

open_cd4naive <- rownames(da_peaks[da_peaks$avg_logFC > 0.5, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_logFC < -0.5, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)

head(closest_genes_cd4naive)

# set plotting order
levels(pbmc) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK CD56bright","NK CD56Dim","pre-B",'pro-B',"pDC","DC","CD14 Mono",'CD16 Mono')

CoveragePlot(
  object = pbmc,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)

