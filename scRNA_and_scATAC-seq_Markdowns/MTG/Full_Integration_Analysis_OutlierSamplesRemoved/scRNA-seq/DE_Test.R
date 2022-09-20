
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
library(sciplot)
library(reshape2)
library(MAST)
library(Matrix.utils)
library(optparse)
library(stats)
library(fgsea)

option_list = list(
  make_option(c("-f", "--SeuratObject"), type="character", default=NULL,
              help="Path to Seurat Object", metavar="character"),
    make_option(c("-c", "--ClusterIdent"), type="character", default="NULL",
              help="Cluster Identity", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

SeuratObjectPath <- opt$SeuratObject

Cluster_Ident <- opt$ClusterIdent

Cluster_Ident

SeuratObject <- readRDS(SeuratObjectPath)

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons", `2` = "Oligodendrocytes", `3` = "GLU_Neurons", `4` = "GLU_Neurons", `5` = "GLU_Neurons", `6` = "Astrocytes", `7` = "GLU_Neurons", `8` = "GLU_Neurons",`9` = "GLU_Neurons",`10` = "GABA_Neurons", `11` = "GLU_Neurons", `12` = "Microglia",`13` = "GABA_Neurons",`14` = "OPCs",`15` = "GABA_Neurons", `16`= "GABA_Neurons", `17`="GABA_Neurons", `18`="Endothelial_Cells", `19`="Endothelial_Cells", `21` = "GABA_Neurons", `22` = "GLU_Neurons", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GLU_Neurons",`27` = "GLU_Neurons", `30` = "Cajal_Retzius", `32` = "GLU_Neurons", `33` = "GLU_Neurons",`34` = "GABA_Neurons", `35` = "GLU_Neurons", `36` = "GABA_Neurons", `38` = "GABA_Neurons", `41` = "GABA_Neurons", `42` = "GABA_Neurons",`43` = "Cajal_Retzius", `44` = "GLU_Neurons", `45` = "GLU_Neurons", `47` = "TEMRA_T_Cells", `48` = "GABA_Neurons", `50` = "GABA_Neurons", `54` = "GLU_Neurons", `66` = "Unknown_Cluster_66")

C2_Gene_Sets <- gmtPathways("~/databases/GSEA_Tables/c2.cp.v7.5.1.symbols.gmt")

dir.create("Files/GSEA")
dir.create("Files/DE_Genes")

ID_Order <- unique(SeuratObject@meta.data$sample_ID)

SeuratObject@meta.data$grouped_case <- SeuratObject@meta.data$case

SeuratObject@meta.data$grouped_case <- gsub("PD|ILB","PD_and_ILB",SeuratObject@meta.data$grouped_case)

TotalGenes <- dim(SeuratObject)[1]

TotalPathways <- length(C2_Gene_Sets)

mkfilter <- function(cmatrixline) {
        sum(cmatrixline > 0)/length(cmatrixline)*100
}

Ident_Object <- subset(SeuratObject, idents = Cluster_Ident)
Ident_Object <- Ident_Object[apply(Ident_Object@assays$RNA@counts,1,mkfilter) >= 20,]
Ident_Object@meta.data$DetRate <- as.vector(scale(colSums(Ident_Object@assays$RNA@counts > 0)))
Idents(Ident_Object) <- "case"
FilteredGeneCount <- dim(Ident_Object)[1]

dir.create("Files/DE_Tests")

metadataDF <- Ident_Object@meta.data

write.table(metadataDF, file = paste("Files/DE_Tests/",Cluster_Ident,"DETests_Metadata_Table.tsv", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")

