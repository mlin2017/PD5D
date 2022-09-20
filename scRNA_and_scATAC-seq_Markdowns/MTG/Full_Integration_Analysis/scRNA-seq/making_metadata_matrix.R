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
  #assign(x = paste(batch,".data", sep = ""), value = Read10X(data.dir = paste("~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/",i,"/scRNA-seq/Files/cellranger_matrices/cellranger_aggr_matrices/filtered_feature_bc_matrix/", sep = "")))
  SeuratObject <- CreateSeuratObject(counts = Read10X(data.dir = paste("~/PD5D_Repo/PD5D/scRNA_and_scATAC-seq_Markdowns/",i,"/scRNA-seq/Files/cellranger_matrices/cellranger_aggr_matrices/filtered_feature_bc_matrix/", sep = "")),
              project = paste(batch,"_data",sep=""),
              min.cells = 3)
  SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^MT-")
  SeuratObject <- subset(SeuratObject, subset = percent.mt < 5)
  assign(paste(batch,"_Object",sep=""),subset(SeuratObject, subset = nFeature_RNA > 200 & nFeature_RNA < 9378.428))
#assign(paste(batch,"_Object",sep=""),subset(SeuratObject, subset = nFeature_RNA > 200 & nFeature_RNA < 9378.428 & percent.mt < 5))
  
}

Batches <- ls()[grep("[[:digit:]]+_Object",ls())]

createmetadata <- function(batch) {
  batchdf <- get(batch)
  tempcellcodes <- as.data.frame(rownames(batchdf@meta.data))
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

write.table(metadata, "Files/metadata_dataframe.tsv", sep = "\t", quote = FALSE)


