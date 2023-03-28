rm(list = ls())
library(Seurat)

#Make object
twentyfour.data <- Read10X(data.dir = "../24M_ham/filtered_feature_bc_matrix")
twentyfour <- CreateSeuratObject(twentyfour.data, project = "twentyfour",min.cells=3,min.features=200)
twentyfour

#Pre-processing 
twentyfour [["percent.mt"]] <- PercentageFeatureSet(twentyfour, pattern = "^MT-")

VlnPlot(twentyfour, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(twentyfour, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(twentyfour, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(twentyfour, feature1 = "nFeature_RNA", feature2 = "percent.mt")

plot1
plot3

#subset cells
twentyfour <- subset(twentyfour, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<20)

#normalize
twentyfour <- NormalizeData(twentyfour, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
twentyfour <- FindVariableFeatures(twentyfour, selection.method = "vst", nfeatures = 2000)
#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(twentyfour), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(twentyfour)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#scaling
all.genes <- rownames(twentyfour)
twentyfour <- ScaleData(twentyfour, features = all.genes)

#linear dimensional reduction
twentyfour <- RunPCA(twentyfour, features = VariableFeatures(object = twentyfour))
print(twentyfour[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(twentyfour, dims = 1:2, reduction = "pca")
DimPlot(twentyfour, reduction = "pca")
DimHeatmap(twentyfour, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(twentyfour, dims = 1:15, cells = 500, balanced = TRUE)


ElbowPlot(twentyfour,ndims = 50)


#easier loop res 3->10
for (i in 3:12){
  #cluster!
  twentyfour <- FindNeighbors(twentyfour, dims = 1:25)
  twentyfour <- FindClusters(twentyfour, resolution = (i/10))
  head(Idents(twentyfour), 5)
  
  #non-linear reduction to umap
  twentyfour <- RunUMAP(twentyfour, dims = 1:25)
  DimPlot(twentyfour, reduction = "umap", label = TRUE)
}

Idents(twentyfour)<-"RNA_snn_res.0.3"
DimPlot(twentyfour, reduction = "umap", label = TRUE, pt.size = 1)

####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(twentyfour, file = "../24M_ham/twentyfour.rds")