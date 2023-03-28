rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)

#Make object
Adipo.data <- Read10X(data.dir = "../Adipo")
Adipo <- CreateSeuratObject(Adipo.data, project = "Adipo",min.cells=3,min.features=200)
Adipo

#Pre-processing 
Adipo [["percent.mt"]] <- PercentageFeatureSet(Adipo, pattern = "^MT-")

VlnPlot(Adipo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Adipo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Adipo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Adipo, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1 + plot2 + plot3
plot1
plot3

#subset cells
Adipo <- subset(Adipo, subset = nFeature_RNA > 2500 & nFeature_RNA < 9500 & percent.mt<20)

#normalize
Adipo <- NormalizeData(Adipo, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
Adipo <- FindVariableFeatures(Adipo, selection.method = "vst", nfeatures = 2000)
#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Adipo), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Adipo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#scaling
all.genes <- rownames(Adipo)
Adipo <- ScaleData(Adipo, features = all.genes)

#linear dimensional reduction
Adipo <- RunPCA(Adipo, features = VariableFeatures(object = Adipo))
print(Adipo[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Adipo, dims = 1:2, reduction = "pca")
DimPlot(Adipo, reduction = "pca")
DimHeatmap(Adipo, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Adipo, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(Adipo,ndims = 50)


#easier loop res 3->10
for (i in 3:12){
  #cluster!
  Adipo <- FindNeighbors(Adipo, dims = 1:30)
  Adipo <- FindClusters(Adipo, resolution = (i/10))
  head(Idents(Adipo), 5)
  
  #non-linear reduction to umap
  Adipo <- RunUMAP(Adipo, dims = 1:30)
  DimPlot(Adipo, reduction = "umap", label = TRUE)
}

Idents(Adipo)<-"RNA_snn_res.0.9"
DimPlot(Adipo, reduction = "umap", label = TRUE, pt.size = 1)


####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(Adipo, file = "../Adipo.rds")