rm(list = ls())
library(Seurat)

#Make object
cuff_72_1.11.22.data <- Read10X(data.dir = "../72_cuff/filtered_feature_bc_matrix")
cuff_72_1.11.22 <- CreateSeuratObject(cuff_72_1.11.22.data, project = "cuff_72_1.11.22",min.cells=3,min.features=200)
cuff_72_1.11.22

#Pre-processing 
cuff_72_1.11.22 [["percent.mt"]] <- PercentageFeatureSet(cuff_72_1.11.22, pattern = "^MT-")

VlnPlot(cuff_72_1.11.22, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(cuff_72_1.11.22, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot3 <- FeatureScatter(cuff_72_1.11.22, feature1 = "nFeature_RNA", feature2 = "percent.mt")

plot1
plot3

#subset cells
cuff_72_1.11.22 <- subset(cuff_72_1.11.22, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt<10)

#normalize
cuff_72_1.11.22 <- NormalizeData(cuff_72_1.11.22, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
cuff_72_1.11.22 <- FindVariableFeatures(cuff_72_1.11.22, selection.method = "vst", nfeatures = 2000)
#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cuff_72_1.11.22), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cuff_72_1.11.22)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#scaling
all.genes <- rownames(cuff_72_1.11.22)
cuff_72_1.11.22 <- ScaleData(cuff_72_1.11.22, features = all.genes)

#linear dimensional reduction
cuff_72_1.11.22 <- RunPCA(cuff_72_1.11.22, features = VariableFeatures(object = cuff_72_1.11.22))
print(cuff_72_1.11.22[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cuff_72_1.11.22, dims = 1:2, reduction = "pca")
DimPlot(cuff_72_1.11.22, reduction = "pca")
DimHeatmap(cuff_72_1.11.22, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(cuff_72_1.11.22, dims = 1:15, cells = 500, balanced = TRUE)


ElbowPlot(cuff_72_1.11.22,ndims = 50)

#easier loop res 3->10
for (i in 3:12){
  #cluster!
  cuff_72_1.11.22 <- FindNeighbors(cuff_72_1.11.22, dims = 1:35)
  cuff_72_1.11.22 <- FindClusters(cuff_72_1.11.22, resolution = (i/10))
  head(Idents(cuff_72_1.11.22), 5)
  
  #non-linear reduction to umap
  cuff_72_1.11.22 <- RunUMAP(cuff_72_1.11.22, dims = 1:35)
  DimPlot(cuff_72_1.11.22, reduction = "umap", label = TRUE)
}

Idents(cuff_72_1.11.22)<-"RNA_snn_res.0.8"
DimPlot(cuff_72_1.11.22, reduction = "umap", label = TRUE, pt.size = 1)

####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(cuff_72_1.11.22, file = "..72_cuff.rds")
