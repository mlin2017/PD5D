library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(dplyr)
library(psych)
library(pheatmap)
library(harmony)
#library(clusterProfiler)
library(DOSE)
library(GOSemSim)
library(enrichplot)
library(stringr)
library(glmpca)
library(SeuratWrappers)
library(Matrix)
library(stats)

for (i in list.files(path = "~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/Midbrain/Multiome/", pattern = "batch[[:digit:]]+_ASAP_Midbrain_Multiome_[[:digit:]]+")){
  aggr <-  read.delim(paste("~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/Midbrain/Multiome/",i,"/snRNA-seq/Files/cellranger_matrices/cellranger_aggr_matrices/aggregation.csv", sep = ""), stringsAsFactors = FALSE, sep = ",")
  all_samples <- aggr$library_id
  batchnumber <- str_extract(i,"batch[[:digit:]]+")
for (id in all_samples){
  filepath <- paste("~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/Midbrain/Multiome/",i,"/snRNA-seq/Files/cellranger_matrices/",id,"/filtered_feature_bc_matrix/",sep = "")
  filetable <- as.data.frame(filepath)
  write.table(filetable, file = "Files/filetable.tsv", sep = "\t")
  SampleMatrix <- Read10X(data.dir = filepath)
  SampleSeurat <- CreateSeuratObject(counts = SampleMatrix,
              project = id,
              min.cells = 3)
  metadata <- as.data.frame(SampleMatrix@Dimnames[[2]])
  colnames(metadata) <- "barcodes"
  metadata$samples <- as.vector(rep(id, ncol(SampleMatrix)))
  metadata$case <- as.vector(rep(aggr[aggr$library_id %in% id,]$case, ncol(SampleMatrix)))
  metadata$sex <- as.vector(rep(aggr[aggr$library_id %in% id,]$sex, ncol(SampleMatrix)))
  metadata$RIN <- as.vector(rep(aggr[aggr$library_id %in% id,]$RIN, ncol(SampleMatrix)))
  metadata$PMI <- as.vector(rep(aggr[aggr$library_id %in% id,]$PMI, ncol(SampleMatrix)))
  metadata$age <- as.vector(rep(aggr[aggr$library_id %in% id,]$age, ncol(SampleMatrix)))
  metadata$batch <- as.vector(rep(batchnumber, ncol(SampleMatrix)))
  SampleSeurat@meta.data$sample_id <- metadata$samples
  SampleSeurat@meta.data$case <- metadata$case
  SampleSeurat@meta.data$batch <- metadata$batch
  SampleSeurat@meta.data$sex <- metadata$sex
  SampleSeurat@meta.data$RIN <- metadata$RIN
  SampleSeurat@meta.data$PMI <- metadata$PMI
  SampleSeurat@meta.data$age <- metadata$age
  age_bracket <- as.vector(cut(SampleSeurat@meta.data$age, c(50,60,70,80,90,100,110)))
  age_bracket <- gsub("\\(|]","",age_bracket)
  age_bracket <- gsub(",","-",age_bracket)
  SampleSeurat@meta.data$age_bracket <- age_bracket
  assign(as.character(as.name(id)), SampleSeurat)
}
}

SampleObjects <- ls()[grep("BN",ls())]

SeuratObject <- merge(get(SampleObjects[1]), 
                 y = lapply(SampleObjects[2:length(SampleObjects)],get), 
                 add.cell.ids = SampleObjects, 
                 project = "batches1-15Integration")

SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^MT-")

nfeature_RNA <- SeuratObject@meta.data$nFeature_RNA
mean(nfeature_RNA)
MAD <- mad(nfeature_RNA, center = median(nfeature_RNA))
threeMAD <- (MAD*3)+median(SeuratObject@meta.data$nFeature_RNA)

SeuratObject <- subset(SeuratObject, subset = nFeature_RNA > 200 & nfeature_RNA < threeMAD & percent.mt < 5)

#saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/Midbrain/Batches1-15_Midbrain_Combined_Matrix_altAB.rds")

cellbarcodesdf <- data.frame(cbind(rownames(SeuratObject@meta.data)[100000:100200],colnames(SeuratObject@assays$RNA@counts)[100000:100200]))

colnames(cellbarcodesdf) <- c("metadata","countmatrix")

write.table(cellbarcodesdf, file = "Files/cell_barcodes_df.txt", sep = "\t", quote = FALSE)

SeuratObject <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(SeuratObject)

SeuratObject <- ScaleData(SeuratObject, features = all.genes, verbose = FALSE)

#finding the top 30 principal components for cells
SeuratObject <- RunGLMPCA(SeuratObject, features=SeuratObject@assays$RNA@var.features, L = 40)

SeuratObject <- RunHarmony(SeuratObject, group.by.vars = c("batch","sex"), plot_convergence = TRUE, reduction = "glmpca", theta = c(1,1))

harmony_embeddings <- Embeddings(SeuratObject, 'harmony')

SeuratObject_ElbowPlot <- ElbowPlot(SeuratObject, reduction = "harmony", ndims = 40)

ggsave2("Figures/SeuratObject_ElbowPlot.pdf", SeuratObject_ElbowPlot, device = "pdf", width = 4, height = 4, units = "in")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/Midbrain/Batches1-15_Midbrain_reducedpcareducedharmony_Part1.rds")
