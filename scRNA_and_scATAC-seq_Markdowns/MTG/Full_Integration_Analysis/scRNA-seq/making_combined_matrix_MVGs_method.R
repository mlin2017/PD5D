#Workflow up to elbow plot to determine cutoff, and saving intermediate seurat
#object as .rds file in /n/scratch3/users/j/jap0606/FullIntegration

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

dir.create("Files")
dir.create("Figures")

for (i in list.files(path = "~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/", pattern = "batch[[:digit:]]+_ASAP_snRNA-seq_[[:digit:]]+")) {
  batch <- str_extract(i,"batch[[:digit:]]+")
  assign(x = paste(batch,"samples", sep = ""), value = read.csv(file.path("~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/",i,"/scRNA-seq/Files/cellranger_matrices/cellranger_aggr_matrices", "aggregation.csv"), stringsAsFactors=F, colClasses = c("sex"="character")))
  assign(x = paste(batch,".data", sep = ""), value = Read10X(data.dir = paste("~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/",i,"/scRNA-seq/Files/cellranger_matrices/cellranger_aggr_matrices/filtered_feature_bc_matrix/", sep = "")))
}

Batches <- ls()[grep("[[:digit:]]+.data",ls())]

createmetadata <- function(batch) {
  batchdf <- get(batch)
  tempcellcodes <- as.data.frame(batchdf@Dimnames[[2]])
  colnames(tempcellcodes) <- "barcodes"
  tempcellcodes$libcodes <- as.factor(gsub(pattern=".+-", replacement="", tempcellcodes$barcodes))
  samptable <- paste("batch",str_extract(batch,"[[:digit:]]+"),"samples",sep="")
  tempcellcodes$samples <- as.vector(get(paste("batch",str_extract(batch,"[[:digit:]]+"),"samples",sep=""))$library_id[tempcellcodes$libcodes])
  tempcellcodes$case <- as.vector(get(paste("batch",str_extract(batch,"[[:digit:]]+"),"samples",sep=""))$case[tempcellcodes$libcodes])
  tempcellcodes$sex <- as.vector(get(paste("batch",str_extract(batch,"[[:digit:]]+"),"samples",sep=""))$sex[tempcellcodes$libcodes])
  tempcellcodes$RIN <- as.vector(get(paste("batch",str_extract(batch,"[[:digit:]]+"),"samples",sep=""))$RIN[tempcellcodes$libcodes])
  tempcellcodes$PMI <- as.vector(get(paste("batch",str_extract(batch,"[[:digit:]]+"),"samples",sep=""))$PMI[tempcellcodes$libcodes])
  tempcellcodes$age <- as.vector(get(paste("batch",str_extract(batch,"[[:digit:]]+"),"samples",sep=""))$age[tempcellcodes$libcodes])
  tempcellcodes$batch <- str_extract(batch,"batch[[:digit:]]")
  return(tempcellcodes)
}

metadata <- bind_rows(lapply(Batches,createmetadata), .id = "column_label")

SeuratObject <- CreateSeuratObject(counts = do.call(cbind,lapply(Batches,get)),
                            project = "FullIntegration_MTG",
                            min.cells = 100, min.features = 200)

SeuratObject@meta.data$sample_id <- metadata$samples
SeuratObject@meta.data$case <- metadata$case
SeuratObject@meta.data$batch <- metadata$batch
SeuratObject@meta.data$sex <- metadata$sex
SeuratObject@meta.data$RIN <- metadata$RIN
SeuratObject@meta.data$PMI <- metadata$PMI
SeuratObject@meta.data$age <- metadata$age
age_bracket <- cut(SeuratObject@meta.data$age, c(50,60,70,80,90,100,110))
age_bracket <- gsub("\\(|]","",age_bracket)
age_bracket <- gsub(",","-",age_bracket)
SeuratObject@meta.data$age_bracket <- age_bracket

SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^MT-")

library(stats)
nfeature_RNA <- SeuratObject@meta.data$nFeature_RNA
mean(nfeature_RNA)
MAD <- mad(nfeature_RNA, center = median(nfeature_RNA))
threeMAD <- (MAD*3)+median(SeuratObject@meta.data$nFeature_RNA)

SeuratObject <- subset(SeuratObject, subset = nFeature_RNA > 200 & nfeature_RNA < threeMAD & percent.mt < 5)

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration/FullIntegration_CombinedFilteredMatrix.rds")


