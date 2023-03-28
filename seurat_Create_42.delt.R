rm(list = ls())
library(Seurat)

#Make object
seur.data <- Read10X(data.dir = "../42_delt/filtered_feature_bc_matrix")
seur <- CreateSeuratObject(seur.data, project = "seur",min.cells=3,min.features=200)
seur

#Pre-processing 
seur [["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-")

VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seur, feature1 = "nFeature_RNA", feature2 = "percent.mt")

plot1
plot3

#subset cells
seur <- subset(seur, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15)

#normalize
seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seur), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seur)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#scaling
all.genes <- rownames(seur)
seur <- ScaleData(seur, features = all.genes)

#linear dimensional reduction
seur <- RunPCA(seur, features = VariableFeatures(object = seur))
print(seur[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seur, dims = 1:2, reduction = "pca")
DimPlot(seur, reduction = "pca")
DimHeatmap(seur, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seur, dims = 1:15, cells = 500, balanced = TRUE)


ElbowPlot(seur,ndims = 50)


#easier loop res 3->10
for (i in 3:12){
  #cluster!
  seur <- FindNeighbors(seur, dims = 1:30)
  seur <- FindClusters(seur, resolution = (i/10))
  head(Idents(seur), 5)
  
  #non-linear reduction to umap
  seur <- RunUMAP(seur, dims = 1:30)
  DimPlot(seur, reduction = "umap", label = TRUE)
}

Idents(seur)<-"RNA_snn_res.0.8"
DimPlot(seur, reduction = "umap", label = TRUE, pt.size = 1)

####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(seur, file = "../42_delt.rds")
