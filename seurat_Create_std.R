rm(list = ls())
library(Seurat)

#Make object
std.data <- Read10X(data.dir = "../Std")
std <- CreateSeuratObject(std.data, project = "std",min.cells=3,min.features=200)
std

#Pre-processing 
std [["percent.mt"]] <- PercentageFeatureSet(std, pattern = "^MT-")

VlnPlot(std, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(std, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(std, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(std, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1 + plot2 + plot3
plot1

#subset cells
std <- subset(std, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt<10)

#normalize
std <- NormalizeData(std, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
std <- FindVariableFeatures(std, selection.method = "vst", nfeatures = 2000)
#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(std), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(std)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#scaling
all.genes <- rownames(std)
std <- ScaleData(std, features = all.genes)

#linear dimensional reduction
std <- RunPCA(std, features = VariableFeatures(object = std))
print(std[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(std, dims = 1:2, reduction = "pca")
DimPlot(std, reduction = "pca")
DimHeatmap(std, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(std, dims = 1:15, cells = 500, balanced = TRUE)


ElbowPlot(std,ndims = 50)

#cluster!
std <- FindNeighbors(std, dims = 1:30)
std <- FindClusters(std, resolution = 0.6)
head(Idents(std), 5)

#non-linear reduction to umap
std <- RunUMAP(std, dims = 1:30)
DimPlot(std, reduction = "umap", label = TRUE, pt.size = 1)

####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(std, file = "../Std.rds")
