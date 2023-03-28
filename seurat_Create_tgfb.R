rm(list = ls())
library(Seurat)

#Make object tgfb = fibro
TGFB.data <- Read10X(data.dir = "../TGFB")
TGFB <- CreateSeuratObject(TGFB.data, project = "TGFB",min.cells=3,min.features=200)
TGFB

#Pre-processing 
TGFB [["percent.mt"]] <- PercentageFeatureSet(TGFB, pattern = "^MT-")

VlnPlot(TGFB, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(TGFB, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TGFB, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(TGFB, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1 + plot2 + plot3
plot1
plot3

#subset cells
TGFB <- subset(TGFB, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt<6)

#normalize
TGFB <- NormalizeData(TGFB, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
TGFB <- FindVariableFeatures(TGFB, selection.method = "vst", nfeatures = 2000)
#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(TGFB), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(TGFB)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#scaling
all.genes <- rownames(TGFB)
TGFB <- ScaleData(TGFB, features = all.genes)

#linear dimensional reduction
TGFB <- RunPCA(TGFB, features = VariableFeatures(object = TGFB))
print(TGFB[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(TGFB, dims = 1:2, reduction = "pca")
DimPlot(TGFB, reduction = "pca")
DimHeatmap(TGFB, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(TGFB, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(TGFB,ndims = 50)


#easier loop res 3->10
for (i in 3:9){
  #cluster!
  TGFB <- FindNeighbors(TGFB, dims = 1:30)
  TGFB <- FindClusters(TGFB, resolution = (i/10))
  head(Idents(TGFB), 5)
  
  #non-linear reduction to umap
  TGFB <- RunUMAP(TGFB, dims = 1:30)
  DimPlot(TGFB, reduction = "umap", label = TRUE)
}

for (res in c("RNA_snn_res.0.3","RNA_snn_res.0.4","RNA_snn_res.0.5","RNA_snn_res.0.6","RNA_snn_res.0.7","RNA_snn_res.0.8","RNA_snn_res.0.9")){
Idents(TGFB) <- res
print(DimPlot(TGFB, reduction = "umap", label = TRUE))
}

Idents(TGFB) <- "RNA_snn_res.0.8"
print(DimPlot(TGFB, reduction = "umap", label = TRUE, pt.size = 1))

####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(TGFB, file = "../TGFB.rds")
