
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

SeuratObject <- readRDS(SeuratObjectPath)

SeuratObject <- RenameIdents(SeuratObject, `0` = "GLU_Neurons", `1` = "GLU_Neurons", `2` = "Oligodendrocytes",
                      `3` = "GLU_Neurons", `4` = "Oligodendrocytes", `5` = "GABA_Neurons",
                      `6` = "GLU_Neurons", `7` = "Astrocytes", `8` = "GLU_Neurons",`9` = "GLU_Neurons",
                      `10` = "GABA_Neurons", `11` = "Microglia",`12` = "GABA_Neurons",
                      `13` = "GLU_Neurons",`14` = "GABA_Neurons",
                      `15` = "OPCs", `16`= "GLU_Neurons", `17`="Endothelial", `18`="Endothelial", `19`="GLU_Neurons",  `20`="GABA_Neurons", `21` = "Oligo_GLU_Neurons", `22` = "GLU_Neurons", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GLU_Neurons",`27` = "GABA_Neurons",`28` = "GLU_Neurons",`29` = "GLU_Neurons",`30` = "Oligo_GABA_Neurons",`31` = "GABA_Neurons",`32` = "GLU_Neurons",`33` = "GABA_Neurons",`34` = "Astrocytes", `35` = "GLU_Neurons")

ID_Order <- unique(SeuratObject@meta.data$sample_id)

#sex <- c("M","M","M","M","M","F","F","M","M","M","F","F","M","M","M","F","F","M","F","M","F","F","M","M","F","M","F","M","F","M","M")

#SeuratObject@meta.data$sex <- SeuratObject@meta.data %>% group_by(sample_id) %>% mutate(Sex = rep(sex[match(unique(sample_id), ID_Order)], length(sample_id))) %>% .$Sex

mkfilter <- function(cmatrixline) {
        sum(cmatrixline > 0)/length(cmatrixline)*100
} 

Ident_Object <- subset(SeuratObject, idents = Cluster_Ident)
Cluster_Ident
Ident_Object <- Ident_Object[apply(Ident_Object@assays$RNA@counts,1,mkfilter) >= 20,]
Ident_Object@meta.data$DetRate <- as.vector(scale(colSums(Ident_Object@assays$RNA@counts)))
Idents(Ident_Object) <- "case"

#PD vs HC
Ident.PDvsHC_Prog.Markers <- FindMarkers(Ident_Object, ident.1 = "PD", ident.2 = "HC", verbose = FALSE, test.use = "MAST", latent.vars = c("sex","DetRate","batch","RIN","PMI","age"))
Ident.PDvsHC_Prog.Markers$gene <- rownames(Ident.PDvsHC_Prog.Markers)
Ident.PDvsHC_Prog.Markers <-
        Ident.PDvsHC_Prog.Markers[Ident.PDvsHC_Prog.Markers$p_val_adj <= 0.05,]

if (length(Ident.PDvsHC_Prog.Markers$gene) != 0) {
        avg.Ident_Object_PDvsHC <- as.data.frame(AverageExpression(Ident_Object, verbose = FALSE)$RNA)
        avg.Ident_Object_PDvsHC$gene <- rownames(avg.Ident_Object_PDvsHC)
        
        avg.Ident_Object_progmarkers_PDvsHC <- avg.Ident_Object_PDvsHC[avg.Ident_Object_PDvsHC$gene %in% unique(Ident.PDvsHC_Prog.Markers$gene),]
        avg.Ident_Object_progmarkers_PDvsHC <- avg.Ident_Object_progmarkers_PDvsHC[match(Ident.PDvsHC_Prog.Markers$gene,avg.Ident_Object_progmarkers_PDvsHC$gene),]
        Ident.PDvsHC_Prog.Markers$PD_mean <- avg.Ident_Object_progmarkers_PDvsHC$PD
        Ident.PDvsHC_Prog.Markers$HC_mean <- avg.Ident_Object_progmarkers_PDvsHC$HC
        Ident.PDvsHC_Prog.Markers$Status <- "Upregulated"
        Ident.PDvsHC_Prog.Markers$Status[Ident.PDvsHC_Prog.Markers$avg_log2FC < 0] <- "Downregulated"
        Ident.PDvsHC_Prog.Markers$Status <- factor(Ident.PDvsHC_Prog.Markers$Status, c("Upregulated","Downregulated"))
        
        write.csv(Ident.PDvsHC_Prog.Markers,file = paste("Files/DE_Genes/All_SigGenes_",Cluster_Ident,"_Markers_HC_vs_PD.csv",sep = ""), quote = FALSE)
}



# ILB vs HC


Ident.ILBvsHC_Prog.Markers <- FindMarkers(Ident_Object, ident.1 = "ILB", ident.2 = "HC", verbose = FALSE, test.use = "MAST", latent.vars = c("sex","DetRate"))
Ident.ILBvsHC_Prog.Markers$gene <- rownames(Ident.ILBvsHC_Prog.Markers)
Ident.ILBvsHC_Prog.Markers <- Ident.ILBvsHC_Prog.Markers[Ident.ILBvsHC_Prog.Markers$p_val_adj <= 0.05,]

if (length(Ident.ILBvsHC_Prog.Markers$gene) != 0){
        avg.Ident_Object_ILBvsHC <- as.data.frame(AverageExpression(Ident_Object, verbose = FALSE)$RNA)
        avg.Ident_Object_ILBvsHC$gene <- rownames(avg.Ident_Object_ILBvsHC)
        
        avg.Ident_Object_progmarkers_ILBvsHC <- avg.Ident_Object_ILBvsHC[avg.Ident_Object_ILBvsHC$gene %in% unique(Ident.ILBvsHC_Prog.Markers$gene),]
        avg.Ident_Object_progmarkers_ILBvsHC <- avg.Ident_Object_progmarkers_ILBvsHC[match(Ident.ILBvsHC_Prog.Markers$gene,avg.Ident_Object_progmarkers_ILBvsHC$gene),]
        Ident.ILBvsHC_Prog.Markers$ILB_mean <- avg.Ident_Object_progmarkers_ILBvsHC$ILB
        Ident.ILBvsHC_Prog.Markers$HC_mean <- avg.Ident_Object_progmarkers_ILBvsHC$HC
        Ident.ILBvsHC_Prog.Markers$Status <- "Upregulated"
        Ident.ILBvsHC_Prog.Markers$Status[Ident.ILBvsHC_Prog.Markers$avg_log2FC < 0] <- "Downregulated"
        Ident.ILBvsHC_Prog.Markers$Status <- factor(Ident.ILBvsHC_Prog.Markers$Status, c("Upregulated","Downregulated"))
        
        write.csv(Ident.ILBvsHC_Prog.Markers,file = paste("Files/DE_Genes/All_SigGenes_",Cluster_Ident,"_Markers_HC_vs_ILB.csv",sep = ""), quote = FALSE)
}

# PD vs ILB

Ident.PDvsILB_Prog.Markers <- FindMarkers(Ident_Object, ident.1 = "PD", ident.2 = "ILB", verbose = FALSE, test.use = "MAST", latent.vars = c("sex","DetRate"))
Ident.PDvsILB_Prog.Markers$gene <- rownames(Ident.PDvsILB_Prog.Markers)
Ident.PDvsILB_Prog.Markers <-
        Ident.PDvsILB_Prog.Markers[Ident.PDvsILB_Prog.Markers$p_val_adj <= 0.05,]

if (length(Ident.PDvsILB_Prog.Markers$gene) != 0) {
        avg.Ident_Object_PDvsILB <- as.data.frame(AverageExpression(Ident_Object, verbose = FALSE)$RNA)
        avg.Ident_Object_PDvsILB$gene <- rownames(avg.Ident_Object_PDvsILB)
        
        avg.Ident_Object_progmarkers_PDvsILB <- avg.Ident_Object_PDvsILB[avg.Ident_Object_PDvsILB$gene %in% unique(Ident.PDvsILB_Prog.Markers$gene),]
        avg.Ident_Object_progmarkers_PDvsILB <- avg.Ident_Object_progmarkers_PDvsILB[match(Ident.PDvsILB_Prog.Markers$gene,avg.Ident_Object_progmarkers_PDvsILB$gene),]
        Ident.PDvsILB_Prog.Markers$PD_mean <- avg.Ident_Object_progmarkers_PDvsILB$PD
        Ident.PDvsILB_Prog.Markers$ILB_mean <- avg.Ident_Object_progmarkers_PDvsILB$ILB
        Ident.PDvsILB_Prog.Markers$Status <- "Upregulated"
        Ident.PDvsILB_Prog.Markers$Status[Ident.PDvsILB_Prog.Markers$avg_log2FC < 0] <- "Downregulated"
        Ident.PDvsILB_Prog.Markers$Status <- factor(Ident.PDvsILB_Prog.Markers$Status, c("Upregulated","Downregulated"))
        
        write.csv(Ident.PDvsILB_Prog.Markers,file = paste("Files/DE_Genes/All_SigGenes_",Cluster_Ident,"_Markers_ILB_vs_PD.csv",sep = ""), quote = FALSE)
}



