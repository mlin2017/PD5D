options(error = function() traceback(2))
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

SeuratObject <- RenameIdents(SeuratObject, `1` = "GLU_Neurons", `2` = "Oligodendrocytes", `3` = "GLU_Neurons", `4` = "GLU_Neurons", `5` = "GLU_Neurons", `6` = "Astrocytes", `7` = "GLU_Neurons", `8` = "GLU_Neurons",`9` = "GLU_Neurons",`10` = "GABA_Neurons", `11` = "GLU_Neurons", `12` = "Microglia",`13` = "GABA_Neurons",`14` = "OPCs",`15` = "GABA_Neurons", `16`= "GABA_Neurons", `17`="GABA_Neurons", `18`="Endothelial_Cells", `19`="Endothelial_Cells", `21` = "GABA_Neurons", `22` = "GLU_Neurons", `23` = "GLU_Neurons", `24` = "GLU_Neurons", `25` = "GLU_Neurons", `26` = "GLU_Neurons",`27` = "GLU_Neurons", `30` = "GABA_Neurons", `32` = "GLU_Neurons", `33` = "GLU_Neurons",`34` = "GABA_Neurons", `35` = "GLU_Neurons", `36` = "GABA_Neurons", `38` = "GABA_Neurons", `41` = "GABA_Neurons", `42` = "GABA_Neurons",`43` = "GABA_Neurons", `44` = "GLU_Neurons", `45` = "GLU_Neurons", `47` = "TEMRA_T_Cells", `48` = "GABA_Neurons", `50` = "GABA_Neurons", `54` = "GLU_Neurons", `66` = "Unknown_Cluster_66")

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

#The filter threshold is a %, DO NOT CHANGE WITHOUT GOOD REASON!
Ident_Object <- subset(SeuratObject, idents = Cluster_Ident)
Ident_Object <- Ident_Object[apply(Ident_Object@assays$RNA@counts,1,mkfilter) >= 20,]
Ident_Object@meta.data$DetRate <- as.vector(scale(colSums(Ident_Object@assays$RNA@counts > 0)))
Idents(Ident_Object) <- "case"
FilteredGeneCount <- dim(Ident_Object)[1]

rm(SeuratObject)

ClusterfGSEA <- function(ClusterIdent, IdentObj, ident1, ident2){
  IdentObj.Markers <- FindMarkers(IdentObj, ident.1 = ident1, ident.2 = ident2, verbose = FALSE, test.use = "MAST", latent.vars = c("batch","RIN","PMI","age","sex","DetRate"), logfc.threshold = 0, min.pct = 0)
  IdentObj.Markers$gene <- rownames(IdentObj.Markers)
  IdentObj.Avg <- as.data.frame(AverageExpression(IdentObj, verbose = FALSE)$RNA)
  IdentObj.Avg$gene <- rownames(IdentObj.Avg)
  ProgMarkersIdentObj.Avg <- IdentObj.Avg[IdentObj.Avg$gene %in% unique(IdentObj.Markers$gene),]
  ProgMarkersIdentObj.Avg <- ProgMarkersIdentObj.Avg[match(IdentObj.Markers$gene,ProgMarkersIdentObj.Avg$gene),]
  generanks <- IdentObj.Markers$avg_log2FC
  names(generanks) <- rownames(IdentObj.Markers)
  generanks <- generanks[order(generanks, decreasing = TRUE)]
  fgseaRes = fgsea(C2_Gene_Sets, stats=generanks, minSize=5, maxSize=Inf, nPermSimple=10000)
  fgseaResAll <- fgseaRes		
  fgseaResAll$leadingEdge = vapply(fgseaResAll$leadingEdge, paste, collapse = ", ", character(1L))		
  write.table(fgseaResAll, file = paste("Files/GSEA/",ClusterIdent,"_",ident1,"_vs_",ident2,"_GSEA_AllGenesets.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  fgseaRes$BF_pval <- 0.05/TotalPathways
  fgseaRes$leadingEdge = vapply(fgseaRes$leadingEdge, paste, collapse = ", ", character(1L))
  fgseaResFilterBH <- fgseaRes[fgseaRes$padj <= 0.05,]
  if (nrow(fgseaResFilterBH) > 0) {
  write.table(fgseaResFilterBH, file = paste("Files/GSEA/",ClusterIdent,"_",ident1,"_vs_",ident2,"_GSEA_SigBHCorrection_Genesets.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  }
  fgseaResFilter <- fgseaRes[fgseaRes$pval <= fgseaRes$BF_pval,]
  if (nrow(fgseaResFilter) > 0) {
  write.table(fgseaResFilter, file = paste("Files/GSEA/",ClusterIdent,"_",ident1,"_vs_",ident2,"_GSEA_SigBFCorrection_Genesets.tsv",sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  }
  IdentObj.Markers$ident1_mean <- ProgMarkersIdentObj.Avg[[ident1]]
  IdentObj.Markers$ident2_mean <- ProgMarkersIdentObj.Avg[[ident2]]
  colnames(IdentObj.Markers)[colnames(IdentObj.Markers) == "ident1_mean"] <- paste(ident1,"_mean",sep="")
  colnames(IdentObj.Markers)[colnames(IdentObj.Markers) == "ident2_mean"] <- paste(ident2,"_mean",sep="")
  IdentObj.Markers$Status <- "Upregulated"
  IdentObj.Markers$Status[IdentObj.Markers$avg_log2FC < 0] <- "Downregulated"
  IdentObj.Markers$Status <- factor(IdentObj.Markers$Status, c("Upregulated","Downregulated"))
  IdentObj.Markers$FilteredGeneCount <- FilteredGeneCount
  IdentObj.Markers$BF_p_val <- 0.05/TotalGenes
  IdentObj.Markers$BH_p_val <- p.adjust(IdentObj.Markers$p_val, method = "BH", n = FilteredGeneCount)
  write.csv(IdentObj.Markers, file = paste("Files/DE_Genes/AllGenes_",ClusterIdent,"_Markers_",ident1,"_vs_",ident2,".csv",sep = ""), quote = FALSE)
}

#PD vs HC
ClusterDE <- function(ClusterIdent, IdentObj, ident1, ident2){
  Ident.PDvsHC_Prog.Markers <- FindMarkers(IdentObj, ident.1 = ident1, ident.2 = ident2, verbose = FALSE, test.use = "MAST", latent.vars = c("batch","RIN","PMI","age","sex","DetRate"), logfc.threshold = 0, min.pct = 0)
  Ident.PDvsHC_Prog.Markers$gene <- rownames(Ident.PDvsHC_Prog.Markers)
  if (length(Ident.PDvsHC_Prog.Markers$gene) != 0) {
        avg.Ident_Object_PDvsHC <- as.data.frame(AverageExpression(Ident_Object, verbose = FALSE)$RNA)
        avg.Ident_Object_PDvsHC$gene <- rownames(avg.Ident_Object_PDvsHC)
        
        avg.Ident_Object_progmarkers_PDvsHC <- avg.Ident_Object_PDvsHC[avg.Ident_Object_PDvsHC$gene %in% unique(Ident.PDvsHC_Prog.Markers$gene),]
        avg.Ident_Object_progmarkers_PDvsHC <- avg.Ident_Object_progmarkers_PDvsHC[match(Ident.PDvsHC_Prog.Markers$gene,avg.Ident_Object_progmarkers_PDvsHC$gene),]
        Ident.PDvsHC_Prog.Markers$ident1_mean <- avg.Ident_Object_progmarkers_PDvsHC[[ident1]]
        colnames(Ident.PDvsHC_Prog.Markers)[colnames(Ident.PDvsHC_Prog.Markers) %in% "ident1_mean"] <- paste(ident1,"_mean", sep = "")
        Ident.PDvsHC_Prog.Markers$ident2_mean <- avg.Ident_Object_progmarkers_PDvsHC[[ident2]]
        colnames(Ident.PDvsHC_Prog.Markers)[colnames(Ident.PDvsHC_Prog.Markers) %in% "ident2_mean"] <- paste(ident2,"_mean", sep = "")
        Ident.PDvsHC_Prog.Markers$Status <- "Upregulated"
        Ident.PDvsHC_Prog.Markers$Status[Ident.PDvsHC_Prog.Markers$avg_log2FC < 0] <- "Downregulated"
        Ident.PDvsHC_Prog.Markers$Status <- factor(Ident.PDvsHC_Prog.Markers$Status, c("Upregulated","Downregulated"))
        Ident.PDvsHC_Prog.Markers$FilteredGeneCount <- FilteredGeneCount
        #Ident.PDvsHC_Prog.Markers$BF_pval <- 0.05/TotalGenes
        #Ident.PDvsHC_Prog.MarkersBF <- Ident.PDvsHC_Prog.Markers[Ident.PDvsHC_Prog.Markers$p_val <= unique(Ident.PDvsHC_Prog.Markers$BF_pval),]
        Ident.PDvsHC_Prog.MarkersBF <- Ident.PDvsHC_Prog.Markers[Ident.PDvsHC_Prog.Markers$p_val_adj <= 0.05,]
        if (length(Ident.PDvsHC_Prog.MarkersBF$gene) != 0) {
            write.csv(Ident.PDvsHC_Prog.MarkersBF,file = paste("Files/DE_Genes/All_SigGenes_",Cluster_Ident,"_Markers_",ident1,"_vs_",ident2,"_BF.csv",sep = ""), quote = FALSE)
        }
        Ident.PDvsHC_Prog.Markers$BH_pval <- p.adjust(Ident.PDvsHC_Prog.Markers$p_val, method = "BH", n = FilteredGeneCount)
        Ident.PDvsHC_Prog.MarkersBH <- Ident.PDvsHC_Prog.Markers[Ident.PDvsHC_Prog.Markers$BH_pval <= 0.05,]
        if (length(Ident.PDvsHC_Prog.MarkersBH$gene) != 0) {
            write.csv(Ident.PDvsHC_Prog.MarkersBH,file = paste("Files/DE_Genes/All_SigGenes_",Cluster_Ident,"_Markers_",ident1,"_vs_",ident2,"_BH.csv",sep = ""), quote = FALSE)
        }
}
}

#PD vs HC

ClusterfGSEA(Cluster_Ident,Ident_Object,"PD","HC")

ClusterDE(Cluster_Ident,Ident_Object,"PD","HC")

# ILB vs HC

ClusterfGSEA(Cluster_Ident,Ident_Object,"ILB","HC")

ClusterDE(Cluster_Ident,Ident_Object,"ILB","HC")

#PD vs ILB

ClusterfGSEA(Cluster_Ident,Ident_Object,"PD","ILB")

ClusterDE(Cluster_Ident,Ident_Object,"PD","ILB")

#GroupedDisease vs HC

Idents(Ident_Object) <- "grouped_case"

ClusterfGSEA(Cluster_Ident,Ident_Object,"PD_and_ILB","HC")

ClusterDE(Cluster_Ident,Ident_Object,"PD_and_ILB","HC")

