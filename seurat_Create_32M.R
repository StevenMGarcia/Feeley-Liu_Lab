rm(list = ls())
library(Seurat)

#Make object
thirtytwo.data <- Read10X(data.dir = "../32M/filtered_feature_bc_matrix")
thirtytwo <- CreateSeuratObject(thirtytwo.data, project = "thirtytwo",min.cells=3,min.features=200)
thirtytwo

#Pre-processing 
thirtytwo [["percent.mt"]] <- PercentageFeatureSet(thirtytwo, pattern = "^MT-")

VlnPlot(thirtytwo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(thirtytwo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(thirtytwo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(thirtytwo, feature1 = "nFeature_RNA", feature2 = "percent.mt")

plot1
plot3

#subset cells
thirtytwo <- subset(thirtytwo, subset = nFeature_RNA > 600 & nFeature_RNA < 7500 & percent.mt<20)

#normalize
thirtytwo <- NormalizeData(thirtytwo, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
thirtytwo <- FindVariableFeatures(thirtytwo, selection.method = "vst", nfeatures = 2000)
#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(thirtytwo), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(thirtytwo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#scaling
all.genes <- rownames(thirtytwo)
thirtytwo <- ScaleData(thirtytwo, features = all.genes)

#linear dimensional reduction
thirtytwo <- RunPCA(thirtytwo, features = VariableFeatures(object = thirtytwo))
print(thirtytwo[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(thirtytwo, dims = 1:2, reduction = "pca")
DimPlot(thirtytwo, reduction = "pca")
DimHeatmap(thirtytwo, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(thirtytwo, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(thirtytwo,ndims = 50)


#easier loop res 3->12
for (i in 3:12){
  #cluster!
  thirtytwo <- FindNeighbors(thirtytwo, dims = 1:25)
  thirtytwo <- FindClusters(thirtytwo, resolution = (i/10))
  head(Idents(thirtytwo), 5)
  
  #non-linear reduction to umap
  thirtytwo <- RunUMAP(thirtytwo, dims = 1:25)
  DimPlot(thirtytwo, reduction = "umap", label = TRUE)
}

Idents(thirtytwo)<-"RNA_snn_res.0.3"
DimPlot(thirtytwo, reduction = "umap", label = TRUE, pt.size = 1)


####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(thirtytwo, file = "../32M/thirtytwo.rds")