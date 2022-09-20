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

for (i in list.files(path = "~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/", pattern = "batch[[:digit:]]+_ASAP_snRNA-seq_[[:digit:]]+")){
  aggr <-  read.delim(paste("~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/",i,"/scRNA-seq/Files/cellranger_matrices/cellranger_aggr_matrices/aggregation.csv", sep = ""), stringsAsFactors = FALSE, sep = ",")
  healthy_controls <- aggr[aggr$case %in% "HC",]$library_id
  batchnumber <- str_extract(i,"batch[[:digit:]]+")
for (id in healthy_controls){
  filepath <- paste("~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/",i,"/scRNA-seq/Files/cellranger_matrices/",id,"/filtered_feature_bc_matrix/",sep = "")
  SampleMatrix <- Read10X(data.dir = filepath)
  SampleSeurat <- CreateSeuratObject(counts = SampleMatrix,
              project = id,
              min.cells = 3)
  metadata <- as.data.frame(SampleMatrix@Dimnames[[2]])
  colnames(metadata) <- "barcodes"
  metadata$samples <- as.vector(rep(id, ncol(SampleMatrix)))
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
  age_bracket <- cut(SampleSeurat@meta.data$age, c(50,60,70,80,90,100,110))
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
                 project = "FullIntegration_MTG")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls/FullIntegrationOSR_HealthyControls_MTG_Combined_Matrix.rds")

cellbarcodesdf <- data.frame(cbind(rownames(SeuratObject@meta.data)[100000:100200],colnames(SeuratObject@assays$RNA@counts)[100000:100200]))

colnames(cellbarcodesdf) <- c("metadata","countmatrix")

write.table(cellbarcodesdf, file = "Files/cell_barcodes_df.txt", sep = "\t", quote = FALSE)

SeuratObject <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(SeuratObject)

SeuratObject <- ScaleData(SeuratObject, features = all.genes, verbose = FALSE)

#finding the top 30 principal components for cells
SeuratObject <- RunGLMPCA(SeuratObject, features=SeuratObject@assays$RNA@var.features, L = 75)

SeuratObject <- RunHarmony(SeuratObject, group.by.vars = c("sample_id","batch","sex","age_bracket"), plot_convergence = TRUE, reduction = "glmpca", theta = c(0.5,0.5,0.5,0.5))

harmony_embeddings <- Embeddings(SeuratObject, 'harmony')

SeuratObject_ElbowPlot <- ElbowPlot(SeuratObject, reduction = "harmony", ndims = 75)

ggsave2("Figures/SeuratObject_ElbowPlot.pdf", SeuratObject_ElbowPlot, device = "pdf", width = 4, height = 4, units = "in")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved_HealthyControls/FullIntegrationOSR_HealthyControls_MTG_Part1.rds")
