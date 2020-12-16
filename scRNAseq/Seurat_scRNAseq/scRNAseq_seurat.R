library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(dplyr)
library(psych)
library(pheatmap)
library(harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(GOSemSim)
library(enrichplot)

# Load the BRI-318, BRI-319, BRI-320 and BRI-321 dataset
BRI318.data <- Read10X(data.dir = "./BRI318")
BRI319.data <- Read10X(data.dir = "./BRI319")
BRI320.data <- Read10X(data.dir = "./BRI320")
BRI321.data <- Read10X(data.dir = "./BRI321")
B0085.data <- Read10X(data.dir = "./B0085")
H0321.data <- Read10X(data.dir = "./H0321")
H1009.data <- Read10X(data.dir = "./H1009")
H1118.data <- Read10X(data.dir = "./H1118")

# Initialize the Seurat object with the raw (non-normalized data)
genes1<-rownames(BRI318.data)
genes2<-rownames(B0085.data)
genes_common<-intersect(genes1,genes2)
BRI318.data<-BRI318.data[rownames(BRI318.data)%in%genes_common,]
BRI319.data<-BRI319.data[rownames(BRI319.data)%in%genes_common,]
BRI320.data<-BRI320.data[rownames(BRI320.data)%in%genes_common,]
BRI321.data<-BRI321.data[rownames(BRI321.data)%in%genes_common,]
B0085.data<-B0085.data[rownames(B0085.data)%in%genes_common,]
H0321.data<-H0321.data[rownames(H0321.data)%in%genes_common,]
H1009.data<-H1009.data[rownames(H1009.data)%in%genes_common,]
H1118.data<-H1118.data[rownames(H1118.data)%in%genes_common,]

dim(BRI318.data)
dim(H0321.data)

AllMB <- CreateSeuratObject(counts = cbind(BRI318.data,
                                           BRI319.data,
                                           BRI320.data,
                                           BRI321.data,
                                           B0085.data,
                                           H0321.data,
                                           H1009.data,
                                           H1118.data),
                            project = "Midbrain",
                            min.cells = 3)

# AllMB_all <- CreateSeuratObject(counts = cbind(BRI318.data,
#                                                BRI319.data,
#                                                BRI320.data,
#                                                BRI321.data,
#                                                B0085.data,
#                                                H0321.data,
#                                                H1009.data,
#                                                H1118.data),
#                                 project = "Midbrain",
#                                 min.cells = 3)

AllMB@meta.data$sample_ID <- c(rep("BRI318", ncol(BRI318.data)),
                               rep("BRI319", ncol(BRI319.data)),
                               rep("BRI320", ncol(BRI320.data)),
                               rep("BRI321", ncol(BRI321.data)),
                               rep("B0085", ncol(B0085.data)),
                               rep("H0321", ncol(H0321.data)),
                               rep("H1009", ncol(H1009.data)),
                               rep("H1118", ncol(H1118.data))
                               )
AllMB@meta.data$case <- c(rep("PD", ncol(BRI318.data)),
                          rep("HC", ncol(BRI319.data)),
                          rep("PD", ncol(BRI320.data)),
                          rep("HC", ncol(BRI321.data)),
                          rep("PD", ncol(B0085.data)),
                          rep("HC", ncol(H0321.data)),
                          rep("HC", ncol(H1009.data)),
                          rep("PD", ncol(H1118.data))
                          )
AllMB@meta.data$experiment <- c(rep("BRI318", ncol(BRI318.data)),
                          rep("BRI319", ncol(BRI319.data)),
                          rep("BRI320", ncol(BRI320.data)),
                          rep("BRI321", ncol(BRI321.data)),
                          rep("B0085", ncol(B0085.data)),
                          rep("H0321", ncol(H0321.data)),
                          rep("H1009", ncol(H1009.data)),
                          rep("H1118", ncol(H1118.data))
)


AllMB@meta.data$BRI318_vs_rest <- c(rep("BRI318", ncol(BRI318.data)),
                                    rep("rest", ncol(BRI319.data)),
                                    rep("rest", ncol(BRI320.data)),
                                    rep("rest", ncol(BRI321.data)),
                                    rep("rest", ncol(B0085.data)),
                                    rep("rest", ncol(H0321.data)),
                                    rep("rest", ncol(H1009.data)),
                                    rep("rest", ncol(H1118.data)))

AllMB@meta.data$BRI319_vs_rest <- c(rep("rest", ncol(BRI318.data)),
                            rep("BRI319", ncol(BRI319.data)),
                            rep("rest", ncol(BRI320.data)),
                            rep("rest", ncol(BRI321.data)),
                            rep("rest", ncol(B0085.data)),
                            rep("rest", ncol(H0321.data)),
                            rep("rest", ncol(H1009.data)),
                            rep("rest", ncol(H1118.data))
)

AllMB@meta.data$BRI320_vs_rest <- c(rep("rest", ncol(BRI318.data)),
                                    rep("rest", ncol(BRI319.data)),
                                    rep("BRI320", ncol(BRI320.data)),
                                    rep("rest", ncol(BRI321.data)),
                                    rep("rest", ncol(B0085.data)),
                                    rep("rest", ncol(H0321.data)),
                                    rep("rest", ncol(H1009.data)),
                                    rep("rest", ncol(H1118.data))
)

AllMB@meta.data$BRI321_vs_rest <- c(rep("rest", ncol(BRI318.data)),
                                    rep("rest", ncol(BRI319.data)),
                                    rep("rest", ncol(BRI320.data)),
                                    rep("BRI321", ncol(BRI321.data)),
                                    rep("rest", ncol(B0085.data)),
                                    rep("rest", ncol(H0321.data)),
                                    rep("rest", ncol(H1009.data)),
                                    rep("rest", ncol(H1118.data))
)

AllMB@meta.data$B0085_vs_rest <- c(rep("rest", ncol(BRI318.data)),
                                    rep("rest", ncol(BRI319.data)),
                                    rep("rest", ncol(BRI320.data)),
                                    rep("rest", ncol(BRI321.data)),
                                    rep("B0085", ncol(B0085.data)),
                                    rep("rest", ncol(H0321.data)),
                                    rep("rest", ncol(H1009.data)),
                                    rep("rest", ncol(H1118.data))
)

AllMB@meta.data$H0321_vs_rest <- c(rep("rest", ncol(BRI318.data)),
                                   rep("rest", ncol(BRI319.data)),
                                   rep("rest", ncol(BRI320.data)),
                                   rep("rest", ncol(BRI321.data)),
                                   rep("rest", ncol(B0085.data)),
                                   rep("H0321", ncol(H0321.data)),
                                   rep("rest", ncol(H1009.data)),
                                   rep("rest", ncol(H1118.data))
)
AllMB@meta.data$H1009_vs_rest <- c(rep("rest", ncol(BRI318.data)),
                                   rep("rest", ncol(BRI319.data)),
                                   rep("rest", ncol(BRI320.data)),
                                   rep("rest", ncol(BRI321.data)),
                                   rep("rest", ncol(B0085.data)),
                                   rep("rest", ncol(H0321.data)),
                                   rep("H1009", ncol(H1009.data)),
                                   rep("rest", ncol(H1118.data))
)

AllMB@meta.data$H1118_vs_rest <- c(rep("rest", ncol(BRI318.data)),
                                   rep("rest", ncol(BRI319.data)),
                                   rep("rest", ncol(BRI320.data)),
                                   rep("rest", ncol(BRI321.data)),
                                   rep("rest", ncol(B0085.data)),
                                   rep("rest", ncol(H0321.data)),
                                   rep("rest", ncol(H1009.data)),
                                   rep("H1118", ncol(H1118.data))
)





# cell counts in different groups before QC
table(AllMB$sample_ID)
## BRI318 BRI319 BRI320 BRI321 
## 5369   1265    922   1925 
table(AllMB$case)
## HC   PD
## 3190 6291

# perform QC and treat data for PCA
AllMB[["percent.mt"]] <- PercentageFeatureSet(AllMB, pattern = "^MT-")
VlnPlot(AllMB, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
AllMB <- subset(AllMB, subset = nFeature_RNA > 200 & percent.mt < 5)
AllMB <- NormalizeData(AllMB, normalization.method = "LogNormalize", scale.factor = 10000)
AllMB <- FindVariableFeatures(AllMB, selection.method = "vst", nfeatures = 2000)
AllMB
# An object of class Seurat 
# 25381 features across 6549 samples within 1 assay 
# Active assay: RNA (25381 features)

# Run the standard workflow for visualization and clustering
all.genes <- rownames(AllMB)
AllMB <- ScaleData(AllMB, features = all.genes, verbose = FALSE)
AllMB <- RunPCA(AllMB, npcs = 30, verbose = FALSE)
AllMB



# An object of class Seurat 
# 25381 features across 6549 samples within 1 assay 
# Active assay: RNA (25381 features)
# 1 dimensional reduction calculated: pca

# cell counts in different groups after QC
table(AllMB$sample_ID)
## BRI318 BRI319 BRI320 BRI321 
## 3612   1045    680   1212 
table(AllMB$case)
## HC   PD
## 2257 4292

# Examine and visualize PCA results a few different ways
VizDimLoadings(AllMB, dims = 1:2, reduction = "pca")
DimPlot(object = AllMB, reduction = "pca", pt.size = .1, group.by = "case")
VlnPlot(object = AllMB, features = "PC_1", group.by = "case",  pt.size = .1)

## FeatureScatter
plot1 <- FeatureScatter(AllMB, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AllMB, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Run Harmony
AllMB <- RunHarmony(AllMB, group.by.vars = "case", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(AllMB, 'harmony')
harmony_embeddings[1:5, 1:5]

DimPlot(object = AllMB, reduction = "harmony", pt.size = .1, group.by = "case")
VlnPlot(object = AllMB, features = "harmony_1", group.by = "case",  pt.size = .1)

AllMB
# An object of class Seurat 
# 25381 features across 6549 samples within 1 assay 
# Active assay: RNA (25381 features)
# 2 dimensional reductions calculated: pca, harmony

# t-SNE and Clustering with harmony
AllMB <- FindNeighbors(AllMB, reduction = "harmony", dims = 1:20)
AllMB <- FindClusters(AllMB, resolution = 0.5)
#AllMB <- RunTSNE(AllMB, reduction = "harmony", dims = 1:20)
AllMB <- RunUMAP(AllMB, reduction = "harmony", dims = 1:20)

plot1 <- DimPlot(AllMB, reduction = "umap", label = TRUE,pt.size = 0.01)
plot1a <- DimPlot(AllMB, reduction = "tsne", label = TRUE,pt.size = 0.01)

plot2 <- DimPlot(AllMB, reduction = "umap", group.by = "case",pt.size = 0.1)
plot3 <- DimPlot(AllMB, reduction = "umap", group.by = "experiment",pt.size = 0.1)

plot4 <- DimPlot(AllMB, reduction = "tsne", group.by = "BRI318_vs_rest")
plot5 <- DimPlot(AllMB, reduction = "tsne", group.by = "BRI319_vs_rest")
plot6 <- DimPlot(AllMB, reduction = "tsne", group.by = "BRI320_vs_rest")
plot7 <- DimPlot(AllMB, reduction = "tsne", group.by = "BRI321_vs_rest")
plot8 <- DimPlot(AllMB, reduction = "tsne", group.by = "B0085_vs_rest")
plot9 <- DimPlot(AllMB, reduction = "tsne", group.by = "H0321_vs_rest")
plot10 <- DimPlot(AllMB, reduction = "tsne", group.by = "H1009_vs_rest")
plot11 <- DimPlot(AllMB, reduction = "tsne", group.by = "H1118_vs_rest")

plot1
plot1a

plot_grid(plot1, plot2,plot3, ncol = 1)
plot_grid(plot2,plot3, ncol = 1)
plot_grid(plot1,plot3, ncol = 1)
plot_grid(plot1,plot2, ncol = 1)
plot_grid(plot1,plot4, ncol = 1)
plot_grid(plot1,plot5, ncol = 1)
plot_grid(plot1,plot6, ncol = 1)
plot_grid(plot1,plot7, ncol = 1)
plot_grid(plot1,plot8, ncol = 1)
plot_grid(plot1,plot9, ncol = 1)
plot_grid(plot1,plot10, ncol = 1)
plot_grid(plot1,plot11, ncol = 1)
plot_grid(plot2,plot3, ncol = 1)
plot_grid(plot2,plot4, ncol = 1)
plot_grid(plot2,plot5, ncol = 1)
plot_grid(plot2,plot6, ncol = 1)
plot_grid(plot2,plot7, ncol = 1)
plot_grid(plot2,plot8, ncol = 1)
plot_grid(plot2,plot9, ncol = 1)
plot_grid(plot2,plot10, ncol = 1)
plot_grid(plot2,plot11, ncol = 1)



DimPlot(AllMB, reduction = "tsne", split.by = "case", label = TRUE, ncol = 1)

# find markers for every cluster compared to all remaining cells, report only the positive ones
AllMB.markers <- FindAllMarkers(AllMB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#saveRDS(AllMB,"AllMB.rds")
AllMB=readRDS("AllMB.rds")
#saveRDS(AllMB.markers,"AllMB.markers.rds")
AllMB.markers=readRDS("AllMB.markers.rds")
write.csv(AllMB.markers,"AllMB_markers.csv")

unique(AllMB.markers$cluster)

AllMB.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(AllMB.markers, file = "AllMarkers.txt", col.names = TRUE, sep = "\t", quote = FALSE)
top10Markers <- AllMB.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(top10Markers, file = "top10Markers.txt", col.names = TRUE, sep = "\t", quote = FALSE)
features <- unique(top10Markers$gene)
DoHeatmap(AllMB, features = features, size = 2, draw.lines = FALSE, angle = 45,
          hjust = 0.2) + theme(axis.text.y = element_text(size = 5))

# # We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types
# FeaturePlot(AllMB, features = c("AQP4", "PLP1", "CX3CR1", "BCAN", "FLT1", "TH", "SLC6A3", "SLC17A6",
#                                          "SLC32A1"))
# FeaturePlot(AllMB, features = c("RBFOX3"))
# VlnPlot(AllMB, features = c("RBFOX3"), pt.size = 0.2, adjust = TRUE)
# FeaturePlot(AllMB, features = c("SLC17A6", "SLC17A7", "TBR1", "BDNF", "CAMK2G", "SLC1A1"))
# VlnPlot(AllMB, features = c("AQP4", "PLP1", "CX3CR1", "BCAN", "FLT1", "CLDN5", "TH", "SLC6A3", "SLC17A6",
#                                      "SLC32A1"), ncol = 2, pt.size = 0.2, adjust = TRUE)
# VlnPlot(AllMB, features = c("RELN", "GAD1", "GAD2", "SLC32A1", "SLC17A6", "SLC17A7", "CAMK2G", "BDNF"),
#         ncol = 2, pt.size = 0.2, adjust = TRUE)
# FeaturePlot(AllMB, features = c("PDGFRA", "CDH19", "LURAP1L-AS1", "PLP1", "OPALIN", "LINC00844", "RBFOX1",
#                                          "KLK6", "GJB1"))
# VlnPlot(AllMB, features = c("PDGFRA", "CDH19", "LURAP1L-AS1", "PLP1", "OPALIN", "LINC00844", "RBFOX1",
#                                      "CDH20", "KLK6", "GJB1"), ncol = 2, pt.size = 0.2, adjust = TRUE)
# FeaturePlot(AllMB, features = c("OLIG2", "SOX10", "SOX6", "PDGFRA", "BCAN", "VCAN"))
# VlnPlot(AllMB, features = c("OLIG2", "SOX10", "SOX6", "PDGFRA", "BCAN", "VCAN"), ncol = 1, pt.size = 0.2,
#         adjust = TRUE)
# FeaturePlot(AllMB, features = c("CALB1", "SLC17A6", "CCK", "CALB2", "EFNB3"))
# FeaturePlot(AllMB, features = c("SOX6", "SLC6A3", "ZDHHC2", "NRIP3", "SLC18A2"))
# VlnPlot(AllMB, features = c("CALB1", "SLC17A6", "CCK", "CALB2", "EFNB3"), ncol = 1, pt.size = 0.2,
#         adjust = TRUE)
# VlnPlot(AllMB, features = c("SOX6", "SLC6A3", "ZDHHC2", "NRIP3", "SLC18A2"), ncol = 1, pt.size = 0.2,
#         adjust = TRUE)
# FeaturePlot(AllMB, features = c("TH", "SLC6A3", "SLC18A2", "DDC"))
# VlnPlot(AllMB, features = c("TH", "SLC6A3", "SLC18A2", "DDC"), ncol = 1, pt.size = 0.2, adjust = TRUE)
# FeaturePlot(AllMB, features = c("FOXA1", "EN1", "EN2", "LMX1B", "PITX3", "NR4A2"))
# VlnPlot(AllMB, features = c("FOXA1", "EN1", "EN2", "LMX1B", "PITX3", "NR4A2"), ncol = 1, pt.size = 0.2,
#         adjust = TRUE)
# FeaturePlot(AllMB, features = c("NDNF", "SOX6", "SNCG", "ALDH1A1", "SNCA", "IGF1"))
# FeaturePlot(AllMB, features = c("SLC32A1", "GAD1", "GAD2", "SLC17A6", "LPL", "OTX2", "ADCYAP1", "GRP", "VIP"))
# VlnPlot(AllMB, features = c("NDNF", "SOX6", "SNCG", "ALDH1A1", "SNCA", "IGF1"), ncol = 1, pt.size = 0.2,
#         adjust = TRUE)
VlnPlot(AllMB, features = c("SLC32A1", "GAD1", "GAD2", "SLC17A6", "LPL", "OTX2", "ADCYAP1", "GRP", "VIP"),
         ncol = 1)

# 
# Neuron  = ENO2, RBFOX3
# Glutamatergic neurons = SLC17A6, SLC17A7
# GABAergic neurons = SLC32A1, GAD1, GAD2
# Dopaminergic neurons = TH, SLC6A3, SCL18A2
# Astrocytes  = AQP4, GFAP
# Oligodendrocytes  =  PLP1, MBP
# OPCs  =  VCAN, BCAN,
# Microglia = CX3CR1, P2RY12
# Endothelial cells = FLT1, CLDN5


FeaturePlot(AllMB, features = c("ENO2", "RBFOX3","SLC17A6", "SLC17A7","AQP4", "GFAP"))
FeaturePlot(AllMB, features = c("SLC17A6", "SLC17A7"))
FeaturePlot(AllMB, features = c("AQP4", "GFAP"))

FeaturePlot(AllMB, features = c("SLC32A1", "GAD1","GAD2","TH", "SLC6A3","SCL8A2"))
FeaturePlot(AllMB, features = c("TH", "SLC6A3","SCL18A2"))

FeaturePlot(AllMB, features = c("PLP1", "MBP","VCAN", "BCAN","CX3CR1", "P2RY12"))
FeaturePlot(AllMB, features = c("VCAN", "BCAN"))
FeaturePlot(AllMB, features = c("CX3CR1", "P2RY12"))
FeaturePlot(AllMB, features = c("PLP1", "MBP","VCAN", "BCAN","FLT1", "CLDN5"))
FeaturePlot(AllMB, features = c("FLT1", "CLDN5"))


pt1=VlnPlot(AllMB, features = c("ENO2"),pt.size = 0)
pt1=VlnPlot(AllMB, features = c("RBFOX3"),pt.size = 0)
pt1=VlnPlot(AllMB, features = c("SLC17A6"),pt.size = 0)
pt1=VlnPlot(AllMB, features = c("SLC17A7"),pt.size = 0)
pt1=VlnPlot(AllMB, features = c("SLC32A1"),pt.size = 0)
VlnPlot(AllMB, features = c("GAD1"),pt.size = 0)
VlnPlot(AllMB, features = c("GAD2"),pt.size = 0)
VlnPlot(AllMB, features = c("TH"),pt.size = 0)
VlnPlot(AllMB, features = c("SLC6A3"),pt.size = 0)
VlnPlot(AllMB, features = c("SLC18A2"),pt.size = 0)
VlnPlot(AllMB, features = c("AQP4"),pt.size = 0)
VlnPlot(AllMB, features = c("GFAP"),pt.size = 0)
VlnPlot(AllMB, features = c("PLP1"),pt.size = 0)
VlnPlot(AllMB, features = c("MBP"),pt.size = 0)
VlnPlot(AllMB, features = c("VCAN"),pt.size = 0)
VlnPlot(AllMB, features = c("BCAN"),pt.size = 0)
VlnPlot(AllMB, features = c("CX3CR1"),pt.size = 0)
VlnPlot(AllMB, features = c("P2RY12"),pt.size = 0)
VlnPlot(AllMB, features = c("FLT1"),pt.size = 0)
VlnPlot(AllMB, features = c("CLDN5"),pt.size = 0)

VlnPlot(AllMB, features = c("ENO2",
                            "TH",
                            "SLC6A3",
                            "SLC18A2",
                            "SLC17A6",
                            "SLC17A7",
                            "SLC32A1",
                            "GAD1",
                            "GAD2",
                            "AQP4",
                            "GFAP",
                            "PLP1",
                            "OLIG1",
                            "VCAN",
                            "CX3CR1",
                            "P2RY12",
                            "FLT1"),
        pt.size =0,ncol = 1)
AverageExpression(AllMB, features = c("ENO2",
                            "TH",
                            "SLC6A3",
                            "SLC18A2",
                            "SLC17A6",
                            "SLC17A7",
                            "SLC32A1",
                            "GAD1",
                            "GAD2",
                            "AQP4",
                            "GFAP",
                            "PLP1",
                            "OLIG1",
                            "VCAN",
                            "CX3CR1",
                            "P2RY12",
                            "FLT1"),use.counts=TRUE)

install.packages("remotes")
remotes::install_github("sbrn3/disscat")

RidgePlot(AllMB, features=c("ENO2"))
RidgePlot(AllMB, features=c("TH"))
RidgePlot(AllMB, features=c("SLC6A3"))
RidgePlot(AllMB, features=c("SLC18A2"))
RidgePlot(AllMB, features=c("SLC17A6"))
RidgePlot(AllMB, features=c("GAD2"))

VlnPlot(AllMB, features = c("OLIG2", "SOX10", "SOX6", "PDGFRA", "BCAN", "VCAN"), ncol = 1, pt.size = 0.2,
        adjust = TRUE)


VlnPlot(AllMB, features = c("OLIG2", "SOX10", "SOX6", "PDGFRA", "BCAN", "VCAN"), ncol = 1, pt.size = 0.2,
        adjust = TRUE)

AverageExpression(AllMB,features = c("OLIG2", "SOX10", "SOX6"))

#########
# Load ggplot2
library(ggplot2)

# create dummy data
data <- data.frame(
  name=letters[1:5],
  value=sample(seq(4,15),5),
  sd=c(1,0.2,3,2,4)
)

# Most basic error bar
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)

##################
###################
# 0,1,2,3,4,5,6,7,8,10,11,14,17
# 9,12,13,15,16
# 5,11,17,2,0,1,3,7,4,6,8,14,10

# Assigning cell type identity to clusters
AllMB <- RenameIdents(AllMB, `0` = "Oligodendrocytes", `1` = "Oligodendrocytes", `2` = "Astrocytes",
                               `3` = "Oligodendrocytes", `4` = "OPCs", `5` = "Glu_GABA neurons",
                               `6` = "Microglia", `7` = "Oligodendrocytes", `8` = "Endothelial cells",`9` = "Glu_neurons",
                                `10` = "Endothelial cells", `11` = "dopamine_neuron_potentially_DA2",`12` = "novel_Cluster12",
                      `13` = "novel_Cluster13",`14` = "novel_Cluster14",
                                `15` = "GABA_neurons", `16`="dopamine_neuron_potentially_DA1")

# AllMB <- RenameIdents(AllMB,
#                       "Oligodendrocytes" = "0",
#                       "Oligodendrocytes"="1",
#                       "Astrocytes"="2",
#                       "Oligodendrocytes"="3",
#                       "OPCs"="4",
#                       "GABAergic neurons"="5",
#                       "Microglia"="6",
#                       "Oligodendrocytes"="7",
#                       "Endothelial cells"="8",
#                       "Endothelial cells"="10",
#                       "Dopaminergic neurons"="11",
#                       "Endothelial cells"="14",
#                       "Dopaminergic neurons"="17")

AllMB <- RenameIdents(AllMB,'Oligodendrocytes'= '0')

# AllMB <- RenameIdents(AllMB,"Endothelial cells 1" = "Endothelial cells" )

# AllMB <- RenameIdents(AllMB, `0` = "Intermediate_Oligo", `1` = "Mature_Oligo 1", `2` = "Astrocytes",
#                       `3` = "Mature_Oligo 2", `4` = "OPCs", `5` = "Mature_Oligo 3",
#                       `6` = "VTA_Glu_GABA Neuron", `7` = "Microglia 1", `8` = "Endothelial cells 1",
#                       `9` = "Glu Neuron", `10` = "SN_DA Neuron", `11` = "Endothelial cells 2",
#                       `12` = "VTA_DA Neuron", `13` = "Microglia 2", `14` = "GABA Neuron")
# 

DimPlot(AllMB, label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(AllMB, label = TRUE, repel = TRUE, pt.size = 1, split.by = "case", label.size = 2,
        ncol = 1) + NoLegend()

DimPlot(AllMB, label = TRUE, repel = TRUE, pt.size = 1, split.by = "experiment", label.size = 2,
        ncol = 1) + NoLegend()



#saveRDS(AllMB,"AllMB.rds")
AllMB=readRDS("AllMB.rds")
saveRDS(AllMB.markers,"AllMB.markers.rds")
AllMB.markers=readRDS("AllMB.markers.rds")

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="plot_1")




# Check the cell counts of each cluster
table(Idents(AllMB))
## Intermediate_Oligo      Mature_Oligo 1          Astrocytes      Mature_Oligo 2                OPCs      Mature_Oligo 3 
##          1437                1224                 870                 830                 556                 438 
## VTA_Glu_GABA Neuron         Microglia 1 Endothelial cells 1          Glu Neuron        SN_DA Neuron Endothelial cells 2 
##          362                 287                 189                 101                  98                  77 
## VTA_DA Neuron         Microglia 2         GABA Neuron 
##          34                  25                  21

# How many cells in each group?
table(AllMB$case)
## HC   PD  
## 2257 4292 

# How many cells in each cluster of each group? 
table(Idents(AllMB), AllMB$case)
##                       HC   PD
## Intermediate_Oligo   436 1001
## Mature_Oligo 1       401  823
## Astrocytes           298  572
## Mature_Oligo 2       233  597
## OPCs                 123  433
## Mature_Oligo 3        92  346
## VTA_Glu_GABA Neuron  211  151
## Microglia 1          116  171
## Endothelial cells 1  131   58
## Glu Neuron            43   58
## SN_DA Neuron          83   15
## Endothelial cells 2   55   22
## VTA_DA Neuron         13   21
## Microglia 2           12   13
## GABA Neuron           10   11

# What proportion of cells are in each cluster of each group?      
prop.table(table(Idents(AllMB), AllMB$case), margin = 2)
##                          HC          PD
## Intermediate_Oligo  0.193176783 0.233224604
## Mature_Oligo 1      0.177669473 0.191752097
## Astrocytes          0.132033673 0.133271202
## Mature_Oligo 2      0.103234382 0.139095993
## OPCs                0.054497120 0.100885368
## Mature_Oligo 3      0.040762074 0.080615098
## VTA_Glu_GABA Neuron 0.093486930 0.035181733
## Microglia 1         0.051395658 0.039841566
## Endothelial cells 1 0.058041648 0.013513514
## Glu Neuron          0.019051839 0.013513514
## SN_DA Neuron        0.036774479 0.003494874
## Endothelial cells 2 0.024368631 0.005125815
## VTA_DA Neuron       0.005759858 0.004892824
## Microglia 2         0.005316792 0.003028891
## GABA Neuron         0.004430660 0.002562908

# AverageExpression
AverageExp <- AverageExpression(AllMB, features = unique(top10Markers$gene), assays = "RNA")
typeof(AverageExp)
head(AverageExp$RNA)
coorda <- corr.test(AverageExp$RNA, AverageExp$RNA, method = "spearman")
pheatmap(coorda$r)

# Identify differential expressed genes across conditions
## DA markers
FeaturePlot(AllMB, features = c("TH", "SLC6A3", "SLC18A2", "DDC"), cols = c("grey", "red"),
            split.by = "case")
plots <- VlnPlot(AllMB, features = c("TH", "SLC6A3", "SLC18A2", "DDC"), split.by = "case", pt.size = 0.2,
                 adjust = TRUE, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

# DEG in each cell type
AllMB$celltype.case <- paste(Idents(AllMB), AllMB$case, sep = "_")
AllMB$celltype <- Idents(AllMB)
Idents(AllMB) <- "celltype.case"
## Astrocytes
Ast.PD_Prog.Markers <- FindMarkers(AllMB, ident.1 = "Astrocytes_PD", ident.2 = "Astrocytes_HC",
                                   verbose = FALSE)
head(Ast.PD_Prog.Markers, n = 10)
Ast.PD_Prog.Markers$gene <- rownames(Ast.PD_Prog.Markers)
Ast.PD_Prog.Markers <- subset(Ast.PD_Prog.Markers, p_val < 0.05)
top10_pos_Ast.PD_Prog.Markers <- Ast.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_Ast.PD_Prog.Markers <- Ast.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)

Idents(AllMB) <- "celltype"
Astrocytes <- subset(AllMB, idents = "Astrocytes")
Idents(Astrocytes) <- "case"
avg.Ast <- log1p(AverageExpression(Astrocytes, verbose = FALSE)$RNA)
avg.Ast$gene <- rownames(avg.Ast)

posfeatures <- unique(top10_pos_Ast.PD_Prog.Markers$gene)
negfeatures <- unique(top10_neg_Ast.PD_Prog.Markers$gene)

features <- c(posfeatures, negfeatures)
p11 <- ggplot(avg.Ast, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("Astrocytes")
p11 <- LabelPoints(plot = p11, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p11 <- LabelPoints(plot = p11, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                   color = "blue") + theme_cowplot(12)
p12 <- DoHeatmap(Astrocytes, features = features, size = 4, angle = 0,
                 hjust = 0.5) + ggtitle("Astrocytes") + theme(axis.text.y = element_text(size = 8))
VlnPlot(Astrocytes, features = features, adjust = TRUE, pt.size = 0)

## Microglia
Microglia <- subset(AllMB, idents = c("Microglia 1", "Microglia 2"))
Idents(Microglia) <- "case"
Mic.PD_Prog.Markers <- FindMarkers(Microglia, ident.1 = "PD", ident.2 = "HC", verbose = FALSE)
head(Mic.PD_Prog.Markers, n = 10)
Mic.PD_Prog.Markers$gene <- rownames(Mic.PD_Prog.Markers)
Mic.PD_Prog.Markers <- subset(Mic.PD_Prog.Markers, p_val < 0.05)
top10_pos_Mic.PD_Prog.Markers <- Mic.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_Mic.PD_Prog.Markers <- Mic.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)
avg.Mic <- log1p(AverageExpression(Microglia, verbose = FALSE)$RNA)
avg.Mic$gene <- rownames(avg.Mic)
posfeatures_Mic <- unique(top10_pos_Mic.PD_Prog.Markers$gene)
negfeatures_Mic <- unique(top10_neg_Mic.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p13 <- ggplot(avg.Mic, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("Microglia")
p13 <- LabelPoints(plot = p13, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p13 <- LabelPoints(plot = p13, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                   color = "blue") + theme_cowplot(12)
p14 <- DoHeatmap(Microglia, features = features, size = 4, angle = 0,
                 hjust = 0.5) + ggtitle("Microglia") + theme(axis.text.y = element_text(size = 8))
VlnPlot(Microglia, features = features, adjust = TRUE, pt.size = 0)
##############
## SN_DA Neuron
SN_DA.PD_Prog.Markers <- FindMarkers(AllMB, ident.1 = "SN_DA Neuron_PD", ident.2 = "SN_DA Neuron_HC", verbose = FALSE)
head(SN_DA.PD_Prog.Markers, n = 10)
SN_DA.PD_Prog.Markers$gene <- rownames(SN_DA.PD_Prog.Markers)
SN_DA.PD_Prog.Markers <- subset(SN_DA.PD_Prog.Markers, p_val < 0.05)
top10_pos_SN_DA.PD_Prog.Markers <- SN_DA.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_SN_DA.PD_Prog.Markers <- SN_DA.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)

theme_set(theme_cowplot())
Idents(AllMB) <- "celltype"
SN_DA.Neuron <- subset(AllMB, idents = "SN_DA Neuron")
Idents(SN_DA.Neuron) <- "case"
avg.SN_DA.Neuron <- log1p(AverageExpression(SN_DA.Neuron, verbose = FALSE)$RNA)
avg.SN_DA.Neuron$gene <- rownames(avg.SN_DA.Neuron)

posfeatures_SN_DA <- unique(top10_pos_SN_DA.PD_Prog.Markers$gene)
negfeatures_SN_DA <- unique(top10_neg_SN_DA.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p1 <- ggplot(avg.SN_DA.Neuron, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("SN_DA Neuron")
p1 <- LabelPoints(plot = p1, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p1 <- LabelPoints(plot = p1, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                  color = "blue") + theme_cowplot(12)
p2 <- DoHeatmap(SN_DA.Neuron, features = features, size = 4, angle = 0, hjust = 0.5) + ggtitle("SN_DA Neuron")
VlnPlot(SN_DA.Neuron, features = features, adjust = TRUE, pt.size = 0)

DGE_DA1 <- FindAllMarkers(SN_DA.Neuron)
avg_DA1 <- AverageExpression(SN_DA.Neuron, slot = "counts")$RNA
write.table(DGE_DA1, file = "DGE_DA1.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(avg_DA1, file = "avg_DA1.txt", col.names = TRUE, sep = "\t", quote = FALSE)

## VTA_DA Neuron
VTA_DA.PD_Prog.Markers <- FindMarkers(AllMB, ident.1 = "VTA_DA Neuron_PD", ident.2 = "VTA_DA Neuron_HC", verbose = FALSE)
head(VTA_DA.PD_Prog.Markers, n = 10)
VTA_DA.PD_Prog.Markers$gene <- rownames(VTA_DA.PD_Prog.Markers)
VTA_DA.PD_Prog.Markers <- subset(VTA_DA.PD_Prog.Markers, p_val < 0.05)
top10_pos_VTA_DA.PD_Prog.Markers <- VTA_DA.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_VTA_DA.PD_Prog.Markers <- VTA_DA.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)

Idents(AllMB) <- "celltype"
VTA_DA.Neuron <- subset(AllMB, idents = "VTA_DA Neuron")
Idents(VTA_DA.Neuron) <- "case"
avg.VTA_DA.Neuron <- log1p(AverageExpression(VTA_DA.Neuron, verbose = FALSE)$RNA)
avg.VTA_DA.Neuron$gene <- rownames(avg.VTA_DA.Neuron)

posfeatures_VTA_DA <- unique(top10_pos_VTA_DA.PD_Prog.Markers$gene)
negfeatures_VTA_DA <- unique(top10_neg_VTA_DA.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p3 <- ggplot(avg.VTA_DA.Neuron, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("VTA_DA Neuron")
p3 <- LabelPoints(plot = p3, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p3 <- LabelPoints(plot = p3, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                  color = "blue") + theme_cowplot(12)
p4 <- DoHeatmap(VTA_DA.Neuron, features = features, size = 4, angle = 0, hjust = 0.5) + ggtitle("VTA_DA Neuron")
VlnPlot(VTA_DA.Neuron, features = features, adjust = TRUE, pt.size = 0)

DGE_DA2 <- FindAllMarkers(VTA_DA.Neuron)
avg_DA2 <- AverageExpression(VTA_DA.Neuron, slot = "counts")$RNA
write.table(DGE_DA2, file = "DGE_DA2.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(avg_DA2, file = "avg_DA2.txt", col.names = TRUE, sep = "\t", quote = FALSE)

## VTA_Glu_GABA Neuron
VTA_Glu_GABA.PD_Prog.Markers <- FindMarkers(AllMB, ident.1 = "VTA_Glu_GABA Neuron_PD", ident.2 = "VTA_Glu_GABA Neuron_HC", verbose = FALSE)
head(VTA_Glu_GABA.PD_Prog.Markers, n = 10)
VTA_Glu_GABA.PD_Prog.Markers$gene <- rownames(VTA_Glu_GABA.PD_Prog.Markers)
VTA_Glu_GABA.PD_Prog.Markers <- subset(VTA_Glu_GABA.PD_Prog.Markers, p_val < 0.05)
top10_pos_VTA_Glu_GABA.PD_Prog.Markers <- VTA_Glu_GABA.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_VTA_Glu_GABA.PD_Prog.Markers <- VTA_Glu_GABA.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)

Idents(AllMB) <- "celltype"
VTA_Glu_GABA.Neuron <- subset(AllMB, idents = "VTA_Glu_GABA Neuron")
Idents(VTA_Glu_GABA.Neuron) <- "case"
avg.VTA_Glu_GABA.Neuron <- log1p(AverageExpression(VTA_Glu_GABA.Neuron, verbose = FALSE)$RNA)
avg.VTA_Glu_GABA.Neuron$gene <- rownames(avg.VTA_Glu_GABA.Neuron)

posfeatures_VTA_Glu_GABA <- unique(top10_pos_VTA_Glu_GABA.PD_Prog.Markers$gene)
negfeatures_VTA_Glu_GABA <- unique(top10_neg_VTA_Glu_GABA.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p5 <- ggplot(avg.VTA_Glu_GABA.Neuron, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("VTA_Glu_GABA Neuron")
p5 <- LabelPoints(plot = p5, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p5 <- LabelPoints(plot = p5, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                  color = "blue") + theme_cowplot(12)
p6 <- DoHeatmap(VTA_Glu_GABA.Neuron, features = features, size = 4, angle = 0, hjust = 0.5) + ggtitle("VTA_Glu_GABA Neuron")
VlnPlot(VTA_Glu_GABA.Neuron, features = features, adjust = TRUE, pt.size = 0)

## Glu Neuron
Glu.PD_Prog.Markers <- FindMarkers(AllMB, ident.1 = "Glu Neuron_PD", ident.2 = "Glu Neuron_HC", verbose = FALSE)
head(Glu.PD_Prog.Markers, n = 10)
Glu.PD_Prog.Markers$gene <- rownames(Glu.PD_Prog.Markers)
Glu.PD_Prog.Markers <- subset(Glu.PD_Prog.Markers, p_val < 0.05)
top10_pos_Glu.PD_Prog.Markers <- Glu.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_Glu.PD_Prog.Markers <- Glu.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)

Idents(AllMB) <- "celltype"
Glu.Neuron <- subset(AllMB, idents = "Glu Neuron")
Idents(Glu.Neuron) <- "case"
avg.Glu.Neuron <- log1p(AverageExpression(Glu.Neuron, verbose = FALSE)$RNA)
avg.Glu.Neuron$gene <- rownames(avg.Glu.Neuron)

posfeatures_Glu <- unique(top10_pos_Glu.PD_Prog.Markers$gene)
negfeatures_Glu <- unique(top10_neg_Glu.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p7 <- ggplot(avg.Glu.Neuron, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("Glu Neuron")
p7 <- LabelPoints(plot = p7, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p7 <- LabelPoints(plot = p7, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                  color = "blue") + theme_cowplot(12)
p8 <- DoHeatmap(Glu.Neuron, features = features, size = 4, angle = 0, hjust = 0.5) + ggtitle("Glu Neuron")
VlnPlot(Glu.Neuron, features = features, adjust = TRUE, pt.size = 0)

## GABA Neuron
GABA.PD_Prog.Markers <- FindMarkers(AllMB, ident.1 = "GABA Neuron_PD", ident.2 = "GABA Neuron_HC", verbose = FALSE)
head(GABA.PD_Prog.Markers, n = 10)
GABA.PD_Prog.Markers$gene <- rownames(GABA.PD_Prog.Markers)
GABA.PD_Prog.Markers <- subset(GABA.PD_Prog.Markers, p_val < 0.05)
top10_pos_GABA.PD_Prog.Markers <- GABA.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_GABA.PD_Prog.Markers <- GABA.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)

Idents(AllMB) <- "celltype"
GABA.Neuron <- subset(AllMB, idents = "GABA Neuron")
Idents(GABA.Neuron) <- "case"
avg.GABA.Neuron <- log1p(AverageExpression(GABA.Neuron, verbose = FALSE)$RNA)
avg.GABA.Neuron$gene <- rownames(avg.GABA.Neuron)

posfeatures_GABA <- unique(top10_pos_GABA.PD_Prog.Markers$gene)
negfeatures_GABA <- unique(top10_neg_GABA.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p9 <- ggplot(avg.GABA.Neuron, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("GABA Neuron")
p9 <- LabelPoints(plot = p9, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p9 <- LabelPoints(plot = p9, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                  color = "blue") + theme_cowplot(12)
p10 <- DoHeatmap(GABA.Neuron, features = features, size = 4, angle = 0, hjust = 0.5) + ggtitle("GABA Neuron")
VlnPlot(GABA.Neuron, features = features, adjust = TRUE, pt.size = 0)

plot_grid(p1, p3, p5, p7, p9, ncol = 2, axis = "b")
plot_grid(p2, p4, p6, p8, p10, ncol = 2)

## Astrocytes
Ast.PD_Prog.Markers <- FindMarkers(AllMB, ident.1 = "Astrocytes_PD", ident.2 = "Astrocytes_HC",
                                   verbose = FALSE)
head(Ast.PD_Prog.Markers, n = 10)
Ast.PD_Prog.Markers$gene <- rownames(Ast.PD_Prog.Markers)
Ast.PD_Prog.Markers <- subset(Ast.PD_Prog.Markers, p_val < 0.05)
top10_pos_Ast.PD_Prog.Markers <- Ast.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_Ast.PD_Prog.Markers <- Ast.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)

Idents(AllMB) <- "celltype"
Astrocytes <- subset(AllMB, idents = "Astrocytes")
Idents(Astrocytes) <- "case"
avg.Ast <- log1p(AverageExpression(Astrocytes, verbose = FALSE)$RNA)
avg.Ast$gene <- rownames(avg.Ast)

posfeatures_Ast <- unique(top10_pos_Ast.PD_Prog.Markers$gene)
negfeatures_Ast <- unique(top10_neg_Ast.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p11 <- ggplot(avg.Ast, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("Astrocytes")
p11 <- LabelPoints(plot = p11, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p11 <- LabelPoints(plot = p11, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                  color = "blue") + theme_cowplot(12)
p12 <- DoHeatmap(Astrocytes, features = features, size = 4, angle = 0,
                 hjust = 0.5) + ggtitle("Astrocytes") + theme(axis.text.y = element_text(size = 8))
VlnPlot(Astrocytes, features = features, adjust = TRUE, pt.size = 0)

## Microglia
Microglia <- subset(AllMB, idents = c("Microglia 1", "Microglia 2"))
Idents(Microglia) <- "case"
Mic.PD_Prog.Markers <- FindMarkers(Microglia, ident.1 = "PD", ident.2 = "HC", verbose = FALSE)
head(Mic.PD_Prog.Markers, n = 10)
Mic.PD_Prog.Markers$gene <- rownames(Mic.PD_Prog.Markers)
Mic.PD_Prog.Markers <- subset(Mic.PD_Prog.Markers, p_val < 0.05)
top10_pos_Mic.PD_Prog.Markers <- Mic.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_Mic.PD_Prog.Markers <- Mic.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)
avg.Mic <- log1p(AverageExpression(Microglia, verbose = FALSE)$RNA)
avg.Mic$gene <- rownames(avg.Mic)
posfeatures_Mic <- unique(top10_pos_Mic.PD_Prog.Markers$gene)
negfeatures_Mic <- unique(top10_neg_Mic.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p13 <- ggplot(avg.Mic, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("Microglia")
p13 <- LabelPoints(plot = p13, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p13 <- LabelPoints(plot = p13, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                   color = "blue") + theme_cowplot(12)
p14 <- DoHeatmap(Microglia, features = features, size = 4, angle = 0,
                 hjust = 0.5) + ggtitle("Microglia") + theme(axis.text.y = element_text(size = 8))
VlnPlot(Microglia, features = features, adjust = TRUE, pt.size = 0)

## Oligodendrocytes
Oligodendrocytes <- subset(AllMB, idents = c("Intermediate_Oligo", "Mature_Oligo 1", "Mature_Oligo 2",
                                                      "Mature_Oligo 3"))
Idents(Oligodendrocytes) <- "case"
Oli.PD_Prog.Markers <- FindMarkers(Oligodendrocytes, ident.1 = "PD", ident.2 = "HC", verbose = FALSE)
head(Oli.PD_Prog.Markers, n = 10)
Oli.PD_Prog.Markers$gene <- rownames(Oli.PD_Prog.Markers)
Oli.PD_Prog.Markers <- subset(Oli.PD_Prog.Markers, p_val < 0.05)
top10_pos_Oli.PD_Prog.Markers <- Oli.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_Oli.PD_Prog.Markers <- Oli.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)
avg.Oli <- log1p(AverageExpression(Oligodendrocytes, verbose = FALSE)$RNA)
avg.Oli$gene <- rownames(avg.Oli)
posfeatures_Oli <- unique(top10_pos_Oli.PD_Prog.Markers$gene)
negfeatures_Oli <- unique(top10_neg_Oli.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p15 <- ggplot(avg.Oli, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("Oligodendrocytes")
p15 <- LabelPoints(plot = p15, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p15 <- LabelPoints(plot = p15, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                   color = "blue") + theme_cowplot(12)
p16 <- DoHeatmap(Oligodendrocytes, features = features, size = 4, angle = 0,
                 hjust = 0.5) + ggtitle("Oligodendrocytes") + theme(axis.text.y = element_text(size = 8))
VlnPlot(Oligodendrocytes, features = features, adjust = TRUE, pt.size = 0)

## OPCs
OPCs <- subset(AllMB, idents = "OPCs")
Idents(OPCs) <- "case"
OPC.PD_Prog.Markers <- FindMarkers(OPCs, ident.1 = "PD", ident.2 = "HC", verbose = FALSE)
head(OPC.PD_Prog.Markers, n = 10)
OPC.PD_Prog.Markers$gene <- rownames(OPC.PD_Prog.Markers)
OPC.PD_Prog.Markers <- subset(OPC.PD_Prog.Markers, p_val < 0.05)
top10_pos_OPC.PD_Prog.Markers <- OPC.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_OPC.PD_Prog.Markers <- OPC.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)
avg.OPC <- log1p(AverageExpression(OPCs, verbose = FALSE)$RNA)
avg.OPC$gene <- rownames(avg.OPC)
posfeatures_OPC <- unique(top10_pos_OPC.PD_Prog.Markers$gene)
negfeatures_OPC <- unique(top10_neg_OPC.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p17 <- ggplot(avg.OPC, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("OPCs")
p17 <- LabelPoints(plot = p17, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p17 <- LabelPoints(plot = p17, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                   color = "blue") + theme_cowplot(12)
p18 <- DoHeatmap(OPCs, features = features, size = 4, angle = 0,
                 hjust = 0.5) + ggtitle("OPCs") + theme(axis.text.y = element_text(size = 8))
VlnPlot(OPCs, features = features, adjust = TRUE, pt.size = 0)

## Endothelial cells
Endothelial_cells <- subset(AllMB, idents = c("Endothelial cells 1", "Endothelial cells 2"))
Idents(Endothelial_cells) <- "case"
End.PD_Prog.Markers <- FindMarkers(Endothelial_cells, ident.1 = "PD", ident.2 = "HC", verbose = FALSE)
head(End.PD_Prog.Markers, n = 10)
End.PD_Prog.Markers$gene <- rownames(End.PD_Prog.Markers)
End.PD_Prog.Markers <- subset(End.PD_Prog.Markers, p_val < 0.05)
top10_pos_End.PD_Prog.Markers <- End.PD_Prog.Markers %>% top_n(n = 10, wt = avg_logFC)
top10_neg_End.PD_Prog.Markers <- End.PD_Prog.Markers %>% top_n(n = 10, wt = -avg_logFC)
avg.End <- log1p(AverageExpression(Endothelial_cells, verbose = FALSE)$RNA)
avg.End$gene <- rownames(avg.End)
posfeatures_End <- unique(top10_pos_End.PD_Prog.Markers$gene)
negfeatures_End <- unique(top10_neg_End.PD_Prog.Markers$gene)
features <- c(posfeatures, negfeatures)
p19 <- ggplot(avg.End, aes(HC, PD)) + geom_point(size = 0.5) + ggtitle("Endothelial cell")
p19 <- LabelPoints(plot = p19, points = posfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3, color = "red")
p19 <- LabelPoints(plot = p19, points = negfeatures, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                   color = "blue") + theme_cowplot(12)
p20 <- DoHeatmap(Endothelial_cells, features = features, size = 4, angle = 0,
                 hjust = 0.5) + ggtitle("Endothelial cell") + theme(axis.text.y = element_text(size = 8))
VlnPlot(Endothelial_cells, features = features, adjust = TRUE, pt.size = 0)

plot_grid(p11, p13, p15, p17, p19, ncol = 2, axis = "b")
plot_grid(p12, p14, p16, p18, p20, ncol = 2)

# comparison between two DA types of neurons
DA_Neuron <- subset(AllMB, idents = c("DA Neuron 1", "DA Neuron 2"))
DA_Neuron <- subset(AllMB_HC, idents = c("DA Neuron 1", "DA Neuron 2"))

DA_DGE <- FindAllMarkers(DA_Neuron)
data_DA_Profile <- as.data.frame(as.matrix(DA_Neuron@assays$RNA@data))
column_name <- as.data.frame(as.matrix(DA_Neuron@active.ident))
column_name <- t(column_name)
data_DA_Profile <- rbind(data_DA_Profile, column_name)
write.table(data_DA_Profile, file = "data_DA_Profile.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(DA_DGE, file = "DA_DGE.txt", col.names = TRUE, sep = "\t", quote = FALSE)

## DA scatter
DA_Neuron <- RenameIdents(DA_Neuron, `DA 1` = "DA_1", `DA 2` = "DA_2")
avg.DA <- log1p(AverageExpression(DA_Neuron, verbose = FALSE)$RNA)
avg.DA$gene <- rownames(avg.DA)
DA_scatter <- ggplot(avg.DA, aes(DA_1, DA_2)) + geom_point(size = 0.5) + ggtitle("DA Neuron")

DA_DGE <- subset(DA_DGE, p_val_adj < 0.05)
DA_DGE <- split(DA_DGE, DA_DGE$cluster)
DA_1_DGE <- DA_DGE$`DA Neuron 1`
DA_2_DGE <- DA_DGE$`DA Neuron 2`
top10_DA_1.Markers <- DA_1_DGE %>% top_n(n = 10, wt = avg_logFC)
top10_DA_2.Markers <- DA_2_DGE %>% top_n(n = 10, wt = avg_logFC)
DA_1_Markers <- unique(top10_DA_1.Markers$gene)
DA_2_Markers <- unique(top10_DA_2.Markers$gene)

DA_scatter <- LabelPoints(plot = DA_scatter, points = DA_1_Markers, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                          color = "red")
DA_scatter <- LabelPoints(plot = DA_scatter, points = DA_2_Markers, repel = TRUE, xnudge = 0, ynudge = 0, size = 3,
                          color = "blue") + theme_cowplot(12)

DA_scatter

## PD susceptibility and progression genes bar plots
PDsusceptibilitygenes <- c("SNCA", "GBA", "PINK1", "ATP13A2", "VPS35", "LRRK2", "PLA2G6")
PDprogressiongenes <- c("GBA", "RIMS2", "IQCJ???SCHIP1", "TMEM108")

data_barplot <- FetchData(AllMB_HC, vars = c("ident", PDsusceptibilitygenes), slot = "counts")
write.table(data_barplot, file = "data_barplot.txt", col.names = TRUE, sep = "\t", quote = FALSE)

data_barplotMarkers <- read.csv(file = "../harmony/PDsusceptibilitygenes.csv")
data_barplotMarkers$ident <- factor(data_barplotMarkers$ident, levels = levels)
data_barplotMarkers$gene <- factor(data_barplotMarkers$gene, levels = PDsusceptibilitygenes)

ggplot(data_barplotMarkers, aes(x = ident, y = mean, fill = ident)) + 
        geom_bar(aes(x = ident, y = mean), stat = "identity", alpha = 1) + 
        geom_errorbar(aes(x = ident, ymin = mean-SE, ymax = mean+SE, colour = ident), width = 0.4, alpha = 0.9, size = 0.5) + 
        facet_grid(rows = vars(gene), scales = "free_y", switch = "y") + 
        theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 45, face = "bold", vjust = 0.5),
              axis.text.y = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
              strip.background = element_blank(), strip.placement = "outside", 
              strip.text.y = element_text(size = 12, angle = 180, face = "bold")) + NoLegend()

data_barplot <- FetchData(AllMB_HC, vars = c("ident", PDprogressiongenes), slot = "counts")
write.table(data_barplot, file = "data_barplot.txt", col.names = TRUE, sep = "\t", quote = FALSE)

data_barplotMarkers <- read.csv(file = "../harmony/PDprogressiongenes.csv")
data_barplotMarkers$ident <- factor(data_barplotMarkers$ident, levels = levels)
data_barplotMarkers$gene <- factor(data_barplotMarkers$gene, levels = PDprogressiongenes)

ggplot(data_barplotMarkers, aes(x = ident, y = mean, fill = ident)) + 
        geom_bar(aes(x = ident, y = mean), stat = "identity", alpha = 1) + 
        geom_errorbar(aes(x = ident, ymin = mean-SE, ymax = mean+SE, colour = ident), width = 0.4, alpha = 0.9, size = 0.5) + 
        facet_grid(rows = vars(gene), scales = "free_y", switch = "y") + 
        theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 45, face = "bold", vjust = 0.5),
              axis.text.y = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
              strip.background = element_blank(), strip.placement = "outside", 
              strip.text.y = element_text(size = 12, angle = 180, face = "bold")) + NoLegend()


# visualization
my_levels <- c("SN_DA Neuron", "VTA_DA Neuron", "VTA_Glu_GABA Neuron", "Glu Neuron", "GABA Neuron",
               "Intermediate_Oligo", "Mature_Oligo 1", "Mature_Oligo 2", "Mature_Oligo 3", "OPCs", "Astrocytes",
               "Microglia 1", "Microglia 2", "Endothelial cells 1", "Endothelial cells 2")
AllMB@active.ident <- factor(x = AllMB@active.ident, levels = my_levels)

featurestoplot <- c(negfeatures_SN_DA, negfeatures_VTA_DA, negfeatures_VTA_Glu_GABA, negfeatures_Glu, negfeatures_GABA,
                    negfeatures_Ast, negfeatures_Mic, negfeatures_End, posfeatures_SN_DA, posfeatures_VTA_DA,
                    posfeatures_VTA_Glu_GABA, posfeatures_Glu, posfeatures_GABA, posfeatures_Ast, posfeatures_Mic,
                    posfeatures_End)

Idents(AllMB) <- "celltype"
AllMB <- RenameIdents(AllMB, `Intermediate_Oligo` = "Oligodendrocytes", `Mature_Oligo 1` = "Oligodendrocytes",
                      `Mature_Oligo 2` = "Oligodendrocytes", `Mature_Oligo 3` = "Oligodendrocytes",
                      `Microglia 1` = "Microglia", `Microglia 2` = "Microglia",
                      `Endothelial cells 1` = "Endothelial cells", `Endothelial cells 2` = "Endothelial cells")
AllMB@meta.data$majorcelltypes <- Idents(AllMB)

levels_majorcelltypes <- c("SN_DA Neuron", "VTA_DA Neuron", "VTA_Glu_GABA Neuron", "Glu Neuron", "GABA Neuron",
                           "Oligodendrocytes", "OPCs", "Astrocytes", "Microglia", "Endothelial cells")
AllMB@active.ident <- factor(x = AllMB@active.ident, levels = levels_majorcelltypes)

exceptOligo <- subset(AllMB, idents = c("SN_DA Neuron", "VTA_DA Neuron", "VTA_Glu_GABA Neuron",
                                        "Glu Neuron", "GABA Neuron", "Astrocytes", "Microglia", "Endothelial cells"))

exceptOligo <- RenameIdents(exceptOligo, `SN_DA Neuron` = "DA_1", `VTA_DA Neuron` = "DA_2",
                            `VTA_Glu_GABA Neuron` = "Glu_GABA", `Glu Neuron` = "Glu", `GABA Neuron` = "GABA",
                            `Astrocytes` = "Ast", `Microglia` = "Mic", `Endothelial cells` = "End")

exceptOligo$majorcelltype.case <- paste(Idents(exceptOligo), exceptOligo$case, sep = " ")
Idents(exceptOligo) <- "majorcelltype.case"

levelsToPlot <- c("DA_1 HC", "DA_2 HC", "Glu_GABA HC", "Glu HC", "GABA HC", "Ast HC", "Mic HC", "End HC",
                  "DA_1 PD", "DA_2 PD", "Glu_GABA PD", "Glu PD", "GABA PD", "Ast PD", "Mic PD", "End PD")
exceptOligo@active.ident <- factor(x = exceptOligo@active.ident, levels = levelsToPlot)

DoHeatmap(exceptOligo, features = featurestoplot, size = 1.2,
          hjust = 0.5) + theme(axis.text.y = element_text(size = 5)) + ggtitle("Top10 down- and up-regulated genes in each PD pathological cell type")

## Top10 Markers of exceptOligo heatmap
Idents(exceptOligo) <- "majorcelltypes"
levels_majorcelltypes <- c("SN_DA Neuron", "VTA_DA Neuron", "VTA_Glu_GABA Neuron", "Glu Neuron", "GABA Neuron",
                           "Astrocytes", "Microglia", "Endothelial cells")
exceptOligo@active.ident <- factor(x = exceptOligo@active.ident, levels = levels_majorcelltypes)

exceptOligo <- RenameIdents(exceptOligo, `SN_DA Neuron` = "DA Neuron 1", `VTA_DA Neuron` = "DA Neuron 2",
                            `VTA_Glu_GABA Neuron` = "Glu_GABA Neuron")

exceptOligo.Markers <- FindAllMarkers(exceptOligo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10exceptOligoMarkers <- exceptOligo.Markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
features <- unique(top10exceptOligoMarkers$gene)

DoHeatmap(exceptOligo, features = features, size = 2, draw.lines = FALSE, angle = 45,
          hjust = 0) + theme(axis.text.y = element_text(size = 5, face = "bold")) + 
        ggtitle("Top10 Marker genes in each cell type")

ggsave(filename = "Top10 Marker genes.tiff", device = "png", dpi = 2000)


## Top 10 marker genes of Midbrain HC samples
top10MBHC.Markers <- MBHC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
features <- unique(top10MBHC.Markers$gene)

DoHeatmap(AllMB_HC, features = features, size = 2, draw.lines = FALSE, angle = 45,
          hjust = 0) + theme(axis.text.y = element_text(size = 5, face = "bold")) + 
        ggtitle("Top10 marker genes for each cell type in the Midbrain HC samples")

ggsave(filename = "Top10 marker genes for each cell type in the Midbrain HC samples.tiff", device = "png", dpi = 2000)


library(data.table)
data_to_write_out <- as.data.frame(as.matrix(exceptOligo@assays$RNA@scale.data))
data_to_write_out <- data_to_write_out[features,]
column_name <- as.data.frame(as.matrix(exceptOligo@active.ident))
colnames(column_name) <- "cell type"
levelsforheatmap <- c("DA Neuron 1", "DA Neuron 2", "Glu_GABA Neuron", "Glu Neuron", "GABA Neuron", "Astrocytes",
                      "Microglia", "Endothelial cells")
column_name$`cell type` <- factor(x = column_name$`cell type`, levels = levelsforheatmap)
column_name <- t(column_name)
heatmapmatrix <- rbind(data_to_write_out, column_name)
write.table(heatmapmatrix, file = "heatmapmatrix.txt", col.names = TRUE, sep = "\t", quote = FALSE)
heatmapmatrix <- read.csv(file = "../harmony/harmony_plots/DGE/heatmapmatrix.csv")
rownames(heatmapmatrix) <- heatmapmatrix$X
heatmapmatrix <- heatmapmatrix[,-1]

pheatmap(heatmapmatrix, fontsize_row = 5, fontsize = 6,
         cellheight = 6, cluster_cols = FALSE, cluster_rows = TRUE, clustering_distance_rows =  "binary",
         legend = TRUE, treeheight_row = 0, treeheight_col = 0, show_colnames = FALSE, annotation = column_name,
         main = "Top 10 genes for each cell type")


# required figures
AllMB <- RenameIdents(AllMB, `SN_DA Neuron` = "DA Neuron 1", `VTA_DA Neuron` = "DA Neuron 2",
                      `VTA_Glu_GABA Neuron` = "Glu_GABA Neuron")
AllMB$majorcelltypes <- Idents(AllMB)
levels <- c("DA Neuron 1", "DA Neuron 2", "Glu_GABA Neuron", "Glu Neuron", "GABA Neuron", "Astrocytes",
            "Oligodendrocytes", "OPCs", "Microglia", "Endothelial cells")
AllMB@active.ident <- factor(x = AllMB@active.ident, levels = levels)

## Marker genes of major cell types vlnplot
f1 <- VlnPlot(AllMB, features = "TH",
              adjust = TRUE, pt.size = 0) + NoLegend() + theme(axis.text.x = element_blank(),
                                                               axis.text.y = element_blank(), 
                                                               axis.title = element_blank(), 
                                                               plot.title = element_text(hjust = -0.08, size = 8, face = "bold.italic", vjust = -12),
                                                               axis.line.y = element_blank(),
                                                               axis.ticks = element_blank(),
                                                               axis.line.x = element_blank(),
                                                               plot.margin = margin(0,0,0.2,1.5,"cm"))

DimPlot(AllMB, reduction = "tsne", pt.size = 1)
DimPlot(AllMB, reduction = "tsne", pt.size = 1, split.by = "case", label = TRUE, ncol = 1)
DimPlot(AllMB_HC, reduction = "tsne", pt.size = 1)

## Major cell types bar graph
VlnPlot(AllMB_HC, features = c("RBFOX3", "ENO2"), ncol = 1)
Majorcelltype_markergenes <- c("ENO2", "TH", "SLC6A3", "SLC18A2", "SLC17A6", "SLC17A7", "SLC32A1", "GAD1", "GAD2", "AQP4",
                               "GFAP", "PLP1", "OLIG1", "VCAN", "CX3CR1", "P2RY12", "FLT1")

data_barplot <- FetchData(AllMB_HC, vars = c("ident", Majorcelltype_markergenes), slot = "counts")
write.table(data_barplot, file = "data_barplot.txt", col.names = TRUE, sep = "\t", quote = FALSE)

data_barplotMarkers <- read.csv(file = "../harmony/data_barplot.csv")
data_barplotMarkers$ident <- factor(data_barplotMarkers$ident, levels = levels)
data_barplotMarkers$gene <- factor(data_barplotMarkers$gene, levels = Majorcelltype_markergenes)

ggplot(data_barplotMarkers, aes(x = ident, y = mean, fill = ident)) + 
        geom_bar(aes(x = ident, y = mean), stat = "identity", alpha = 1) + 
        geom_errorbar(aes(x = ident, ymin = mean-SE, ymax = mean+SE, colour = ident), width = 0.4, alpha = 0.9, size = 0.5) + 
        facet_grid(rows = vars(gene), scales = "free_y", switch = "y") + 
        theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 45, face = "bold", vjust = 0.5),
              axis.text.y = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
              strip.background = element_blank(), strip.placement = "outside", 
              strip.text.y = element_text(size = 12, angle = 180, face = "bold")) + NoLegend()

## PD related genes expression
Idents(AllMB) <- "case"
AllMB_HC <- subset(AllMB, idents = "HC")
Idents(AllMB_HC) <- "majorcelltypes"

table(AllMB_HC$sample_ID)
table(Idents(AllMB_HC))
table(Idents(AllMB_HC), AllMB_HC$sample_ID)
prop.table(table(Idents(AllMB_HC), AllMB_HC$sample_ID), margin = 2)
table(AllMB_HC$celltype)
table(Idents(AllMB_HC), AllMB_HC$sample_ID)

AllMB_PD <- subset(AllMB, idents = "PD")
Idents(AllMB_PD) <- "majorcelltypes"

VlnPlot(AllMB_HC, features = c("HLA-C",	"HLA-DMA", "HLA-DQA1", "HLA-DRB5", "HLA-DRA", "HLA-DOB", "HLA-DQA2",
                               "HLA-DQB1", "HLA-DRB1", "LRRK2"), ncol = 2, pt.size = 0.1)

VlnPlot(AllMB_HC, features = c("MAPT",	"SNCA",	"VPS35", "UCHL1", "PINK1", "ATP13A2", "HTRA2", "PLA2G6", "RIMS2",
                               "IQCJ-SCHIP1", "TMEM108"), ncol = 3, pt.size = 0.02)

VlnPlot(AllMB_PD, features = c("HLA-C", "HLA-DMA", "HLA-DQA1", "HLA-DRB5", "HLA-DRA", "HLA-DOB", "HLA-DQA2",
                               "HLA-DQB1", "HLA-DRB1", "LRRK2"), ncol = 2, pt.size = 0.1)

VlnPlot(AllMB_PD, features = c("MAPT", "SNCA", "VPS35", "UCHL1", "PINK1", "ATP13A2", "HTRA2", "PLA2G6", "RIMS2",
                               "IQCJ-SCHIP1", "TMEM108"), ncol = 3, pt.size = 0.02)

## heatmap for Midbrain Neuron types
MBHC_Neuron <- subset(AllMB_HC, idents = c("DA Neuron 1", "DA Neuron 2", "Glu_GABA Neuron", "Glu Neuron",
                                           "GABA Neuron"))
AllMB_HC@active.ident <- factor(AllMB_HC@active.ident, levels = levels)
AllMB_HC$majorcelltypes <- Idents(AllMB_HC)

MBHC.markers <- FindAllMarkers(AllMB_HC, only.pos = TRUE, logfc.threshold = 0.25)
write.table(MBHC.markers, file = "MBHCMarkers.txt", col.names = TRUE, sep = "\t", quote = FALSE)
top20MBHC.Markers <- MBHC.markers %>% group_by(cluster) %>% top_n(n = 20, wt = -p_val_adj)
write.table(top20MBHC.Markers, file = "top20MBHC.Markers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

top20MBHC_Neuron_Markers <- top20MBHC.Markers[1:100,]
top20MBHC_Neuron_genes <- unique(top20MBHC_Neuron_Markers$gene)

DAgenes <- c("SOX6", "SNCG", "NDNF", "IGF1", "FOXA2", "LMX1A", "ALDH1A1", "SLC32A1", "SATB1", "CLSTN2", "ADCYAP1",
             "LPL", "OTX2", "VIP", "CHRNA4", "GSG1L", "SNCA", "NTF3", "TH", "SLC6A3", "SLC18A2", "DDC", "SLC17A6",
             "SLC17A7", "GAD1", "GAD2", "ZDHHC2", "NRIP3", "CALB1", "CALB2", "CCK", "EFNB3", "ZFHX3", "TACR3", "NR4A2",
             "GRP")

unique(DAgenes)
heatmapgenes <- unique(c(top20MBHC_Neuron_genes, DAgenes))

Neuronsdataforheatmap <- as.data.frame(as.matrix(MBHC_Neuron@assays$RNA@scale.data))
Neuronsdataforheatmap <- Neuronsdataforheatmap[heatmapgenes,]

column_name <- as.data.frame(as.matrix(MBHC_Neuron@active.ident))
colnames(column_name) <- "cell type"
column_name <- t(column_name)
levelsforheatmap <- c("DA Neuron 1", "DA Neuron 2", "Glu_GABA Neuron", "Glu Neuron", "GABA Neuron")
column_name$`cell type` <- factor(x = column_name$`cell type`, levels = levelsforheatmap)

DAheatmapmatrix <- rbind(Neuronsdataforheatmap, column_name)
write.table(DAheatmapmatrix, file = "DAheatmapmatrix.txt", col.names = TRUE, sep = "\t", quote = FALSE)
DAheatmapmatrix <- read.csv(file = "../harmony/harmony_plots/DGE/DAheatmapmatrix.csv")
rownames(DAheatmapmatrix) <- DAheatmapmatrix$X
DAheatmapmatrix <- DAheatmapmatrix[,-1]

pheatmap(DAheatmapmatrix, fontsize_row = 5, treeheight_col = 0, cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = TRUE, clustering_distance_rows =  "manhattan", annotation_col = column_name,
         annotation_names_col = TRUE, main = "Marker genes of Midbrain Neuronal subtypes")



# enrichment analysis
top10_pos_SN_DA.PD_Prog.Markers$celltype <- rep("DA Neuron 1", nrow(top10_pos_SN_DA.PD_Prog.Markers))
top10_pos_VTA_DA.PD_Prog.Markers$celltype <- rep("DA Neuron 2", nrow(top10_pos_VTA_DA.PD_Prog.Markers))
top10_pos_VTA_Glu_GABA.PD_Prog.Markers$celltype <- rep("Glu_GABA Neuron", nrow(top10_pos_VTA_Glu_GABA.PD_Prog.Markers))
top10_pos_Glu.PD_Prog.Markers$celltype <- rep("Glu Neuron", nrow(top10_pos_Glu.PD_Prog.Markers))
top10_pos_GABA.PD_Prog.Markers$celltype <- rep("GABA Neuron", nrow(top10_pos_GABA.PD_Prog.Markers))

top10_pos_Oli.PD_Prog.Markers$celltype <- rep("Oligodendrocytes", nrow(top10_pos_Oli.PD_Prog.Markers))
top10_pos_OPC.PD_Prog.Markers$celltype <- rep("OPCs", nrow(top10_pos_OPC.PD_Prog.Markers))
top10_pos_Ast.PD_Prog.Markers$celltype <- rep("Astrocytes", nrow(top10_pos_Ast.PD_Prog.Markers))
top10_pos_Mic.PD_Prog.Markers$celltype <- rep("Microglia", nrow(top10_pos_Mic.PD_Prog.Markers))
top10_pos_End.PD_Prog.Markers$celltype <- rep("Endothelial cells", nrow(top10_pos_End.PD_Prog.Markers))

top10_neg_SN_DA.PD_Prog.Markers$celltype <- rep("DA Neuron 1", nrow(top10_neg_SN_DA.PD_Prog.Markers))
top10_neg_VTA_DA.PD_Prog.Markers$celltype <- rep("DA Neuron 2", nrow(top10_neg_VTA_DA.PD_Prog.Markers))
top10_neg_VTA_Glu_GABA.PD_Prog.Markers$celltype <- rep("Glu_GABA Neuron", nrow(top10_neg_VTA_Glu_GABA.PD_Prog.Markers))
top10_neg_Glu.PD_Prog.Markers$celltype <- rep("Glu Neuron", nrow(top10_neg_Glu.PD_Prog.Markers))
top10_neg_GABA.PD_Prog.Markers$celltype <- rep("GABA Neuron", nrow(top10_neg_GABA.PD_Prog.Markers))

top10_neg_Oli.PD_Prog.Markers$celltype <- rep("Oligodendrocytes", nrow(top10_neg_Oli.PD_Prog.Markers))
top10_neg_OPC.PD_Prog.Markers$celltype <- rep("OPCs", nrow(top10_neg_OPC.PD_Prog.Markers))
top10_neg_Ast.PD_Prog.Markers$celltype <- rep("Astrocytes", nrow(top10_neg_Ast.PD_Prog.Markers))
top10_neg_Mic.PD_Prog.Markers$celltype <- rep("Microglia", nrow(top10_neg_Mic.PD_Prog.Markers))
top10_neg_End.PD_Prog.Markers$celltype <- rep("Endothelial cells", nrow(top10_neg_End.PD_Prog.Markers))

Alltop10up_regulatedgenes <- rbind(top10_pos_SN_DA.PD_Prog.Markers, top10_pos_VTA_DA.PD_Prog.Markers,
                                   top10_pos_VTA_Glu_GABA.PD_Prog.Markers, top10_pos_Glu.PD_Prog.Markers,
                                   top10_pos_GABA.PD_Prog.Markers, top10_pos_Oli.PD_Prog.Markers,
                                   top10_pos_OPC.PD_Prog.Markers, top10_pos_Ast.PD_Prog.Markers,
                                   top10_pos_Mic.PD_Prog.Markers, top10_pos_End.PD_Prog.Markers)

Alltop10down_regulatedgenes <- rbind(top10_neg_SN_DA.PD_Prog.Markers, top10_neg_VTA_DA.PD_Prog.Markers,
                                   top10_neg_VTA_Glu_GABA.PD_Prog.Markers, top10_neg_Glu.PD_Prog.Markers,
                                   top10_neg_GABA.PD_Prog.Markers, top10_neg_Oli.PD_Prog.Markers,
                                   top10_neg_OPC.PD_Prog.Markers, top10_neg_Ast.PD_Prog.Markers,
                                   top10_neg_Mic.PD_Prog.Markers, top10_neg_End.PD_Prog.Markers)

Alltop10up_regulatedgenes$celltype <- as.factor(Alltop10up_regulatedgenes$celltype)
Alltop10down_regulatedgenes$celltype <- as.factor(Alltop10down_regulatedgenes$celltype)

## DA Neuron 1 analysis
DA_1_genes <- unique(DA_1_DGE$gene)
DA_1_genes.df <- bitr(DA_1_genes, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db)
DA_1_genes <- unique(DA_1_genes.df$ENTREZID)

DA_1_go <- groupGO(gene = DA_1_pos_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", level = 3, readable = TRUE)
DA_1_go <- groupGO(gene = DA_1_pos_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", level = 4, readable = TRUE)
DA_1_go <- groupGO(gene = DA_1_pos_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", level = 3, readable = TRUE)
DA_1_go <- groupGO(gene = DA_1_pos_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", level = 4, readable = TRUE)

write.table(DA_1_DGE, file = "DA_1_DGE.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(DA_1_go@result, file = "DA_1_go@result.txt", col.names = TRUE, sep = "\t", quote = FALSE)

genelist <- DA_1_DGE[unique(DA_1_genes.df$SYMBOL),2]
names(genelist) <- DA_1_genes
genelist = sort(genelist, decreasing = TRUE)


DA_1_ego <- enrichGO(gene = DA_1_pos_genes, OrgDb = org.Hs.eg.db, ont = "MF", qvalueCutoff = 0.05, readable = TRUE)
DA_1_ego <- enrichGO(gene = DA_1_pos_genes, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
DA_1_ego <- simplify(DA_1_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
DA_1_ego <- enrichGO(gene = DA_1_genes, OrgDb = org.Hs.eg.db, ont = "CC", qvalueCutoff = 0.05, readable = TRUE)

write.table(DA_1_ego@result, file = "DA_1_ego@result.txt", col.names = TRUE, sep = "\t", quote = FALSE)

DA_1_kk <- enrichKEGG(gene = DA_1_pos_genes, organism = "hsa")
DA_1_kk <- setReadable(DA_1_kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.table(DA_1_kk@result, file = "DA_1_kk@result.txt", col.names = TRUE, sep = "\t", quote = FALSE)

f1 <- dotplot(DA_1_ego, title = "upregulated genes GO_ORA")
f3 <- dotplot(DA_1_kk, title = "upregulated genes KEGG_ORA")

## DA Neuron 2 analysis
DA_2_genes <- unique(DA_2_DGE$gene)
DA_2_genes.df <- bitr(DA_2_genes, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db)
DA_2_genes <- unique(DA_2_genes.df$ENTREZID)

DA_2_go <- groupGO(gene = DA_1_neg_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", level = 3, readable = TRUE)
DA_2_go <- groupGO(gene = DA_2_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", level = 4, readable = TRUE)
DA_2_go <- groupGO(gene = DA_2_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", level = 3, readable = TRUE)
DA_2_go <- groupGO(gene = DA_1_neg_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", level = 4, readable = TRUE)

write.table(DA_2_DGE, file = "DA_2_DGE.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(DA_2_go@result, file = "DA_2_go@result.txt", col.names = TRUE, sep = "\t", quote = FALSE)

genelist <- DA_2_DGE[unique(DA_2_genes.df$SYMBOL),2]
names(genelist) <- DA_2_genes
genelist = sort(genelist, decreasing = TRUE)

DA_2_ego <- enrichGO(gene = DA_1_neg_genes, OrgDb = org.Hs.eg.db, ont = "MF", readable = TRUE)
DA_2_ego <- enrichGO(gene = DA_1_neg_genes, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
DA_2_ego <- simplify(DA_2_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
DA_2_ego <- enrichGO(gene = DA_1_neg_genes, OrgDb = org.Hs.eg.db, ont = "CC", readable = TRUE)

DA_2_ego2 <- gseGO(geneList = genelist, OrgDb = org.Hs.eg.db, ont = "MF")

write.table(DA_2_ego@result, file = "DA_2_ego@result.txt", col.names = TRUE, sep = "\t", quote = FALSE)

DA_2_kk <- enrichKEGG(gene = DA_1_neg_genes, organism = "hsa")
DA_2_kk <- setReadable(DA_2_kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.table(DA_2_kk@result, file = "DA_2_kk@result.txt", col.names = TRUE, sep = "\t", quote = FALSE)
 
f2 <- dotplot(DA_2_ego, title = "downregulated genes GO_ORA")
f4 <- dotplot(DA_2_kk, showCategory = 8, title = "downregulated genes KEGG_ORA")
cnetplot(DA_2_kk, foldChange = genelist, circular = TRUE, colorEdge = TRUE)

###
DA_1_DGE <- subset(DA_1_DGE, p_val_adj < 0.05)

DA_1_df <- data.frame(Entrez = names(genelist), FC = genelist)
DA_1_df$group <- "upregulated"
DA_1_df$group[DA_1_df$FC < 0] <- "downregulated"

DA_1_df <- split(DA_1_df, DA_1_df$group)
DA_1_pos <- DA_1_df$upregulated
DA_1_neg <- DA_1_df$downregulated
DA_1_pos_genes <- as.character(DA_1_pos$Entrez)
DA_1_neg_genes <- as.character(DA_1_neg$Entrez)

formula_res <- compareCluster(Entrez~group, data = DA_1_df, fun = "enrichDO")
head(as.data.frame(formula_res))

dotplot(formula_res, showCategory = 30)

plot_grid(f1, f2)

## DA Neuron 2 analysis*
DA_2_allgenes <- subset(MBHC.markers, cluster == "DA Neuron 2")
DA_2_allgenes <- subset(DA_2_allgenes, p_val_adj < 0.05)
DA2allmarkers <- as.character(unique(DA_2_allgenes$gene))
DA2allmarkers.df <- bitr(DA2allmarkers, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db)
DA2allmarkers <- unique(DA2allmarkers.df$ENTREZID)

DA_2AM_ego <- enrichGO(gene = DA2allmarkers, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
DA_2AM_ego <- simplify(DA_2AM_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

DA_2AM_kk <- enrichKEGG(gene = DA2allmarkers, organism = "hsa")
DA_2AM_kk <- setReadable(DA_2AM_kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

DA_2AM_DO <- enrichDO(gene = DA2allmarkers, readable = TRUE)

rownames(DA_2_allgenes) <- DA_2_allgenes$gene
genelist <- DA_2_allgenes[unique(DA2allmarkers.df$SYMBOL),2]
names(genelist) <- DA2allmarkers
genelist = sort(genelist, decreasing = TRUE)

DA_2AM_DO2 <- gseDO(geneList = genelist)

f5 <- dotplot(DA_2AM_DO, showCategory = 20, title = "DA Neuron 2 Allmarkers enrichDO")

## DA Neuron 1 analysis*
DA_1_allgenes <- subset(MBHC.markers, cluster == "DA Neuron 1")
DA_1_allgenes <- subset(DA_1_allgenes, p_val_adj < 0.05)
DA1allmarkers <- as.character(unique(DA_1_allgenes$gene))
DA1allmarkers.df <- bitr(DA1allmarkers, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db)
DA1allmarkers <- unique(DA1allmarkers.df$ENTREZID)

DA_1AM_ego <- enrichGO(gene = DA1allmarkers, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
DA_1AM_ego <- simplify(DA_1AM_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

DA_1AM_kk <- enrichKEGG(gene = DA1allmarkers, organism = "hsa")
DA_1AM_kk <- setReadable(DA_1AM_kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

DA_1AM_DO <- enrichDO(gene = DA1allmarkers, readable = TRUE)

f6 <- dotplot(DA_1AM_DO, showCategory = 20, title = "DA Neuron 1 Allmarkers enrichDO")

plot_grid(f6, f5)


## heatmap for key midbrain neuron genes

keymidbraingenes <- c("SLC18A2", "SLC6A3", "SYT1", "SLC17A6", "NRXN3", "CADM2", "NRXN1", "PARK2", "SNCA", "ALDH1A1",
                      "HSP90AB1", "CKB", "HSP90AA1", "ATP1B1", "MALAT1", "SLC17A7", "GAD2", "SLC32A1", "GAD1")

DoHeatmap(MBHC_Neuron, features = keymidbraingenes, size = 4, draw.lines = FALSE, angle = 45, hjust = 0) + 
        theme(axis.text.y = element_text(size = 10, face = "bold")) + 
        ggtitle("Midbrain Neuronal Subtypes")

Neuronsdataforheatmap <- as.data.frame(as.matrix(MBHC_Neuron@assays$RNA@scale.data))
Neuronsdataforheatmap <- Neuronsdataforheatmap[keymidbraingenes,]

column_name <- as.data.frame(as.matrix(MBHC_Neuron@active.ident))
colnames(column_name) <- "cell type"
levelsforheatmap <- c("DA Neuron 1", "DA Neuron 2", "Glu_GABA Neuron", "GABA Neuron", "Glu Neuron")
column_name$`cell type` <- factor(x = column_name$`cell type`, levels = levelsforheatmap)
column_name <- t(column_name)
column_name <- as.data.frame(column_name)

DAheatmapmatrix <- rbind(Neuronsdataforheatmap, column_name)
write.table(DAheatmapmatrix, file = "DAheatmapmatrix.txt", col.names = TRUE, sep = "\t", quote = FALSE)
DAheatmapmatrix <- read.csv(file = "DAheatmapmatrix.csv" )
rownames(DAheatmapmatrix) <- DAheatmapmatrix$X
DAheatmapmatrix <- DAheatmapmatrix[,-1]

pheatmap(DAheatmapmatrix, fontsize_row = 5, treeheight_col = 0, cluster_cols = FALSE, show_colnames = FALSE,
         cluster_rows = TRUE, clustering_distance_rows =  "euclidean", main = "Midbrain Neuronal Subtypes")

FeaturePlot(MBHC_Neuron, features = "SNCA", label = TRUE)
DimPlot(MBHC_Neuron, reduction = "tsne")

##
MBHC_Neuron.markers <- FindAllMarkers(MBHC_Neuron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10MBHC_Neuron.markers <- MBHC_Neuron.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

levelsforheatmap <- c("DA Neuron 1", "DA Neuron 2", "Glu_GABA Neuron", "Glu Neuron", "GABA Neuron")
MBHC_Neuron@active.ident <- factor(x = MBHC_Neuron@active.ident, levels = levelsforheatmap)

keymidbraingenes <- c("CHRNB3", "LMX1B", "TTC6", "AC007091.1", "TH", "SLC18A2", "SLC6A3", "ALDH1A1", "CANX", "ATP1B1",
                      "SERINC1", "HSP90AA1", "SNCA", "PARK2", "ASIC2", "ROBO1", "TENM2", "SLC17A6", "SLC17A7", "FAT2",
                      "RP11-649A16.1", "MSRA", "FSTL5", "GAD1", "GAD2", "SLC32A1", "CASR", "SHISA8", "DMBX1")

DoHeatmap(MBHC_Neuron, features = keymidbraingenes, size = 4, draw.lines = FALSE, angle = 45, hjust = 0) + 
        theme(axis.text.y = element_text(size = 10, face = "bold")) + 
        ggtitle("Midbrain Neuronal Subtypes")
        

keyDAgenes <- c("MALAT1", "SNCA", "PARK2", "SLC6A3", "ALDH1A1", "HSP90AB1", "HSP90AA1", "CKB", "ATP1B1", "GAD1", "SLC17A7")     
        
DoHeatmap(DA_Neuron, features = keyDAgenes, size = 4, draw.lines = FALSE, angle = 45, hjust = 0) + 
        theme(axis.text.y = element_text(size = 10, face = "bold")) + 
        ggtitle("DA Neuronal Subtypes")        

DA_1_Markers <- FindMarkers(DA_Neuron, ident.1 = "DA Neuron 1", ident.2 = "DA Neuron 2")        
write.table(DA_1_Markers, file = "DA_1_Markers.txt", col.names = TRUE, sep = "\t", quote = FALSE)
DA_2_Markers <- FindMarkers(DA_Neuron, ident.1 = "DA Neuron 2", ident.2 = "DA Neuron 1")            
write.table(DA_2_Markers, file = "DA_2_Markers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

keyDAgenes <- c("PCDH15", "PDE4D", "CSMD1", "MEG3", "TH", "SLC18A2", "SLC6A3", "ALDH1A1", "ATP1B1", "SERINC1",
                "HSP90AA1", "SNCA", "PARK2", "GAD1", "SLC17A7")     


## Adrenoreceptors expression in Major cell types bar graph
VlnPlot(AllMB_HC, features = c("ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2"), ncol = 2)
Adrenoreceptors_genes <- c("ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2")

data_barplot <- FetchData(AllMB_HC, vars = c("ident", Adrenoreceptors_genes), slot = "counts")
write.table(data_barplot, file = "data_barplot.txt", col.names = TRUE, sep = "\t", quote = FALSE)

data_barplotMarkers <- read.csv(file = "Adrenoreceptors_genes HC.csv")
data_barplotMarkers$ident <- factor(data_barplotMarkers$ident, levels = levels)
data_barplotMarkers$gene <- factor(data_barplotMarkers$gene, levels = Adrenoreceptors_genes)

ggplot(data_barplotMarkers, aes(x = ident, y = mean, fill = ident)) + 
        geom_bar(aes(x = ident, y = mean), stat = "identity", alpha = 1) + 
        geom_errorbar(aes(x = ident, ymin = mean-SE, ymax = mean+SE, colour = ident), width = 0.4, alpha = 0.9, size = 0.5) + 
        facet_grid(rows = vars(gene), scales = "free_y", switch = "y") + 
        theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 45, face = "bold", vjust = 0.5),
              axis.text.y = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
              strip.background = element_blank(), strip.placement = "outside", 
              strip.text.y = element_text(size = 12, angle = 180, face = "bold")) + NoLegend() + 
        ggtitle("Adrenoreceptors genes expression in Midbrain HC samples")

        
VlnPlot(AllMB_PD, features = c("ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2"), ncol = 2)

data_barplot <- FetchData(AllMB_PD, vars = c("ident", Adrenoreceptors_genes), slot = "counts")
write.table(data_barplot, file = "data_barplot.txt", col.names = TRUE, sep = "\t", quote = FALSE)

data_barplotMarkers <- read.csv(file = "Adrenoreceptors_genes PD.csv")
data_barplotMarkers$ident <- factor(data_barplotMarkers$ident, levels = levels)
data_barplotMarkers$gene <- factor(data_barplotMarkers$gene, levels = Adrenoreceptors_genes)

ggplot(data_barplotMarkers, aes(x = ident, y = mean, fill = ident)) + 
        geom_bar(aes(x = ident, y = mean), stat = "identity", alpha = 1) + 
        geom_errorbar(aes(x = ident, ymin = mean-SE, ymax = mean+SE, colour = ident), width = 0.4, alpha = 0.9, size = 0.5) + 
        facet_grid(rows = vars(gene), scales = "free_y", switch = "y") + 
        theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 45, face = "bold", vjust = 0.5),
              axis.text.y = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
              strip.background = element_blank(), strip.placement = "outside", 
              strip.text.y = element_text(size = 12, angle = 180, face = "bold")) + NoLegend() + 
        ggtitle("Adrenoreceptors genes expression in Midbrain PD samples")


FeaturePlot(AllMB, features = c("SLC6A1", "SLC6A11"), split.by = "case", label = TRUE)
FeaturePlot(AllMB, features = c("SLC1A3", "SLC1A2", "SLC1A1", "SLC1A6"), split.by = "case", label =TRUE)

VlnPlot(AllMB, features = c("SLC6A1", "SLC6A11"), split.by = "case", ncol = 1)
VlnPlot(AllMB, features = c("SLC1A3", "SLC1A2", "SLC1A1", "SLC1A6"), split.by = "case", ncol = 2)
