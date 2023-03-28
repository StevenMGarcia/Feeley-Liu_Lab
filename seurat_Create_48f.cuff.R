rm(list = ls())
library(Seurat)

#Make object
cuff_48f_2.10.22.data <- Read10X(data.dir = "../48f_cuff/filtered_feature_bc_matrix")
cuff_48f_2.10.22 <- CreateSeuratObject(cuff_48f_2.10.22.data, project = "cuff_48f_2.10.22",min.cells=3,min.features=200)
cuff_48f_2.10.22

#Pre-processing 
cuff_48f_2.10.22 [["percent.mt"]] <- PercentageFeatureSet(cuff_48f_2.10.22, pattern = "^MT-")

#VlnPlot(cuff_48f_2.10.22, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot1 <- FeatureScatter(cuff_48f_2.10.22, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(cuff_48f_2.10.22, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(cuff_48f_2.10.22, feature1 = "nFeature_RNA", feature2 = "percent.mt")

#plot1
plot3

#subset cells
cuff_48f_2.10.22 <- subset(cuff_48f_2.10.22, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15)

#normalize
cuff_48f_2.10.22 <- NormalizeData(cuff_48f_2.10.22, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
cuff_48f_2.10.22 <- FindVariableFeatures(cuff_48f_2.10.22, selection.method = "vst", nfeatures = 2000)

#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cuff_48f_2.10.22), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cuff_48f_2.10.22)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1
plot2

#scaling
all.genes <- rownames(cuff_48f_2.10.22)
cuff_48f_2.10.22 <- ScaleData(cuff_48f_2.10.22, features = all.genes)

#linear dimensional reduction
cuff_48f_2.10.22 <- RunPCA(cuff_48f_2.10.22, features = VariableFeatures(object = cuff_48f_2.10.22))
print(cuff_48f_2.10.22[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cuff_48f_2.10.22, dims = 1:2, reduction = "pca")
DimPlot(cuff_48f_2.10.22, reduction = "pca")
DimHeatmap(cuff_48f_2.10.22, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(cuff_48f_2.10.22,ndims = 50)

#easier loop res 3->12
for (i in 3:12){
  #cluster!
  cuff_48f_2.10.22 <- FindNeighbors(cuff_48f_2.10.22, dims = 1:30)
  cuff_48f_2.10.22 <- FindClusters(cuff_48f_2.10.22, resolution = (i/10))
  head(Idents(cuff_48f_2.10.22), 5)
  
  #non-linear reduction to umap
  cuff_48f_2.10.22 <- RunUMAP(cuff_48f_2.10.22, dims = 1:30)
  DimPlot(cuff_48f_2.10.22, reduction = "umap", label = TRUE)
}

Idents(cuff_48f_2.10.22)<-"RNA_snn_res.0.3"
Idents(cuff_48f_2.10.22)<-"percent.mt"
DimPlot(cuff_48f_2.10.22, reduction = "umap", label = FALSE, pt.size = 1, raster = TRUE)


####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(cuff_48f_2.10.22, file = "../48f_cuff.rds")
