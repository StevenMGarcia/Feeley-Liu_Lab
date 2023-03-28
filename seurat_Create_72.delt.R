rm(list = ls())
library(Seurat)

#Make object
delt_72_1.11.22.data <- Read10X(data.dir = "../72_delt/filtered_feature_bc_matrix")
delt_72_1.11.22 <- CreateSeuratObject(delt_72_1.11.22.data, project = "delt_72_1.11.22",min.cells=3,min.features=200)
delt_72_1.11.22

#Pre-processing 
delt_72_1.11.22 [["percent.mt"]] <- PercentageFeatureSet(delt_72_1.11.22, pattern = "^MT-")

VlnPlot(delt_72_1.11.22, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(delt_72_1.11.22, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(delt_72_1.11.22, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(delt_72_1.11.22, feature1 = "nFeature_RNA", feature2 = "percent.mt")

plot1
plot3

#subset cells
delt_72_1.11.22 <- subset(delt_72_1.11.22, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<10)

#normalize
delt_72_1.11.22 <- NormalizeData(delt_72_1.11.22, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
delt_72_1.11.22 <- FindVariableFeatures(delt_72_1.11.22, selection.method = "vst", nfeatures = 2000)
#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(delt_72_1.11.22), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(delt_72_1.11.22)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#scaling
all.genes <- rownames(delt_72_1.11.22)
delt_72_1.11.22 <- ScaleData(delt_72_1.11.22, features = all.genes)

#linear dimensional reduction
delt_72_1.11.22 <- RunPCA(delt_72_1.11.22, features = VariableFeatures(object = delt_72_1.11.22))
print(delt_72_1.11.22[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(delt_72_1.11.22, dims = 1:2, reduction = "pca")
DimPlot(delt_72_1.11.22, reduction = "pca")
DimHeatmap(delt_72_1.11.22, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(delt_72_1.11.22, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(delt_72_1.11.22,ndims = 50)


#easier loop res 3->10
for (i in 3:12){
  #cluster!
  delt_72_1.11.22 <- FindNeighbors(delt_72_1.11.22, dims = 1:30)
  delt_72_1.11.22 <- FindClusters(delt_72_1.11.22, resolution = (i/10))
  head(Idents(delt_72_1.11.22), 5)
  
  #non-linear reduction to umap
  delt_72_1.11.22 <- RunUMAP(delt_72_1.11.22, dims = 1:30)
  DimPlot(delt_72_1.11.22, reduction = "umap", label = TRUE)
}

Idents(delt_72_1.11.22)<-"RNA_snn_res.0.7"
DimPlot(delt_72_1.11.22, reduction = "umap", label = TRUE, pt.size = 1)


####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(delt_72_1.11.22, file = "../72_delt.rds")
