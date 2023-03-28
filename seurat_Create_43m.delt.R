rm(list = ls())
library(Seurat)

#Make object
delt_43m_1.25.22.data <- Read10X(data.dir = "../43m_delt/filtered_feature_bc_matrix")
delt_43m_1.25.22 <- CreateSeuratObject(delt_43m_1.25.22.data, project = "delt_43m_1.25.22",min.cells=3,min.features=200)
delt_43m_1.25.22

#Pre-processing 
delt_43m_1.25.22 [["percent.mt"]] <- PercentageFeatureSet(delt_43m_1.25.22, pattern = "^MT-")

VlnPlot(delt_43m_1.25.22, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(delt_43m_1.25.22, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(delt_43m_1.25.22, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(delt_43m_1.25.22, feature1 = "nFeature_RNA", feature2 = "percent.mt")

plot1
plot3

#subset cells
delt_43m_1.25.22 <- subset(delt_43m_1.25.22, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15)

#normalize
delt_43m_1.25.22 <- NormalizeData(delt_43m_1.25.22, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
delt_43m_1.25.22 <- FindVariableFeatures(delt_43m_1.25.22, selection.method = "vst", nfeatures = 2000)

#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(delt_43m_1.25.22), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(delt_43m_1.25.22)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1
plot2

#scaling
all.genes <- rownames(delt_43m_1.25.22)
delt_43m_1.25.22 <- ScaleData(delt_43m_1.25.22, features = all.genes)

#linear dimensional reduction
delt_43m_1.25.22 <- RunPCA(delt_43m_1.25.22, features = VariableFeatures(object = delt_43m_1.25.22))
print(delt_43m_1.25.22[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(delt_43m_1.25.22, dims = 1:2, reduction = "pca")
DimPlot(delt_43m_1.25.22, reduction = "pca")
DimHeatmap(delt_43m_1.25.22, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(delt_43m_1.25.22,ndims = 50)

#easier loop res 3->10
for (i in 3:12){
  #cluster!
  delt_43m_1.25.22 <- FindNeighbors(delt_43m_1.25.22, dims = 1:30)
  delt_43m_1.25.22 <- FindClusters(delt_43m_1.25.22, resolution = (i/10))
  head(Idents(delt_43m_1.25.22), 5)
  
  #non-linear reduction to umap
  delt_43m_1.25.22 <- RunUMAP(delt_43m_1.25.22, dims = 1:30)
  DimPlot(delt_43m_1.25.22, reduction = "umap", label = TRUE)
}

Idents(delt_43m_1.25.22)<-"RNA_snn_res.0.7"
DimPlot(delt_43m_1.25.22, reduction = "umap", label = TRUE, pt.size = 1)


####
#SAVE IT!
#SAVE IT!
#SAVE IT!
saveRDS(delt_43m_1.25.22, file = "../43m_delt.rds")