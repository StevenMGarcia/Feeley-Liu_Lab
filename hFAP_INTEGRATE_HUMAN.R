suppressMessages(
  {
    library(dplyr)
    library(Seurat)
    library(patchwork)
    library(tibble)
    library(DoubletFinder)
    library(scCustomize)
    library(monocle3)
    library(SeuratWrappers)
    library(ggplot2)
    library(SeuratDisk)
    library(xlsx)
    library(loupeR)
    library(ggVennDiagram)
    library(stringr)
    library(Seurat)
    library(hdf5r)
    library(loomR)
    library(EnhancedVolcano)
    library(SignacX)
    library(glmGamPoi)
    library(enrichR)
  })
data_path <- "~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Data/"

#mouse_data <- readRDS(paste0(data_path,"GSE200487/GSE200487_regen_data.rds"))
#mouse <- CreateSeuratObject(mouse_data$SCT@data, meta.data = mouse_data$SCT$varMetadata)

# Load Path: Make sure the paths are right for your objets-----------------
cellranger_path <-paste0(data_path,"GSE200487/")
cellranger_data <- paste0(cellranger_path,"cellranger/")
cellranger_folder <- list.dirs(cellranger_data, full.names = FALSE, recursive = FALSE)
#If you want to run a specific subset of sc, then specify a pattern, if not leave as null
obj_path <- paste0(cellranger_path,"Objects/")
proj_name <- "Label_Transfer"

# Runs all folders and outputs seurat objects -----------------------------
for (i in 1:length(cellranger_folder)){
  cellranger <- Read10X(data.dir = paste0(cellranger_data,cellranger_folder[i]))
  assign(cellranger_folder[i],
         obj <- CreateSeuratObject(counts = cellranger, project = cellranger_folder[i], min.cells = 3, min.features = 200))
}


for (i in 1:length(cellranger_folder)){
  
  cellranger <- Read10X(data.dir = paste0(cellranger_data,cellranger_folder[i]))
  # Initialize the Seurat object with the raw (non-normalized data).
  obj <- CreateSeuratObject(counts = cellranger, project = cellranger_folder[i], min.cells = 3, min.features = 200)
  
  # QC ----------------------------------------------------------------------
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
  feature_RNA <- 0.83*max(obj@meta.data$nFeature_RNA)
  #obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < feature_RNA & percent.mt < 10)
  
  # Normalize Data ----------------------------------------------------------
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identify Variable Features ----------------------------------------------
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(obj), 10)  #identify the 10 most highly variable genes
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
  # Scale Data --------------------------------------------------------------
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  
  # Linear Dimensional Reduction --------------------------------------------
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 30)
  
  #Figure out how many dims to include
  
  ElbowPlot(obj,ndims = 30)
  
  # For loop for resolutions 0.3 - 1.2  ------------------------------------------
  for (j in 10:12){
    #cluster!
    obj <- FindNeighbors(obj, dims = 1:30)
    obj <- FindClusters(obj, resolution = (j/10))
    head(Idents(obj), 5)
    
    #non-linear reduction to umap
    obj <- RunUMAP(obj, dims = 1:30)
    DimPlot(obj, reduction = "umap", label = TRUE)
  }
  
  Idents(obj)<-"RNA_snn_res.0.8"
  DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1) #preliminary clustering
  
  # Save merged object ------------------------------------------------------
  saveRDS(obj, file = paste0(obj_path,cellranger_folder[i],".rds"))
  
}


# Load --------------------------------------------------------------------
data_path <- "~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Data/"
cellranger_path <-paste0(data_path,"GSE200487/")
cellranger_data <- paste0(cellranger_path,"cellranger/")
cellranger_folder <- list.dirs(cellranger_data, full.names = FALSE, recursive = FALSE)
obj_path <- paste0(cellranger_path,"Objects/")

obj_list <- list.files(obj_path)
obj_folder <- list.dirs(cellranger_data, full.names = FALSE, recursive = FALSE)
for (i in 1:length(obj_list)){
  assign(cellranger_folder[i],
         readRDS(paste0(obj_path,obj_list[i])))
}

GSE200487 <- merge(RF1,y = c(RF2,GM1,GM2),add.cell.ids = c("RF1", "RF2","GM1", "GM2"), project = "label_transfer")
GSE200487@meta.data$dataset <- "GSE200487"



obj <- GSE200487

#obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
#VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

#obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data$orig.ident)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
#obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  verbose = FALSE, k.weight = 64,
  orig.reduction = "pca", new.reduction = "integrated.cca"
)

obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:40)
obj <- FindClusters(obj, resolution = 1, cluster.name = "cca_clusters", reduction = "integrated.cca")

obj <- RunUMAP(obj, dims = 1:40, reduction = "integrated.cca")
DimPlot(
  obj,
  combine = FALSE, label.size = 2, group.by = "orig.ident", reduction = "umap"
)

DimPlot(obj,label = T)

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

#GSE_markers <- FindAllMarkers(obj)

FeaturePlot_scCustom(obj, features = c("PDGFRA","CD34"))
rm(obj)


GSE_faps <- subset(obj, idents = c("6","10","12"))


saveRDS(GSE_faps, "~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE200487.rds")


saveRDS(GSE_faps, "~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE200487_FAPs.rds")


# Reload ------------------------------------------------------------------
steven_fap <- readRDS("~/Desktop/UCSF/Seurat_Objects/SMG_FAP.rds")
GSE_faps <- readRDS("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE200487_FAPs.rds")
#steven_fap@meta.data$dataset <- "steven"
#DimPlot(steven_fap)
#steven_fap@meta.data$adjusted_celltype <- steven_fap@meta.data$celltype

#levels(steven_fap@meta.data$adjusted_celltype) <- c("CCN3","LAMA2","CD74","DLK1","ATF3","CD55","GLI1","PIEZO2","MYL9","TNMD","TNNT1")
Idents(steven_fap) <- "adjusted_celltype"
DimPlot(steven_fap)
#saveRDS(steven_fap, "~/Desktop/UCSF/Seurat_Objects/SMG_FAP.rds")

merged <- merge(steven_fap,y = GSE_faps, add.cell.ids = c("steven","GSE200487"))

merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$dataset)

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)

merged <- IntegrateLayers(object = merged, method = CCAIntegration, orig.reduction = "pca",new.reduction = "integrated.cca", verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:30)
merged <- FindClusters(merged, resolution = 2, cluster.name = "cca_clusters", reduction = "integrated.cca")

merged <- RunUMAP(merged, reduction = "integrated.cca", dims = 1:30)
DimPlot(merged, group.by = c("dataset", "adjusted_celltype", "cluster"))

DimPlot(merged)

saveRDS(merged, "~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE200487_fap_steven.rds")
merged <- readRDS("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE200487_fap_steven.rds")


Idents(merged) <- "cluster"
Idents(merged) <- "dataset"
steven.ref <- merged
GSE.query <- merged
rm(merged)

DefaultAssay(GSE.query) <- "RNA"
DefaultAssay(steven.ref) <- "RNA"


GSE.query <- NormalizeData(GSE.query)
cell.anchors <- FindTransferAnchors(reference = steven.ref, query = GSE.query, dims = 1:32, reference.reduction = "pca")

predictions <- TransferData(anchorset = cell.anchors, refdata = steven.ref$adjusted_celltype, dims = 1:32)
GSE.query <- AddMetaData(GSE.query, metadata = predictions)


GSE.query$prediction.match <- GSE.query$predicted.id == GSE.query$adjusted_celltype
table(GSE.query$prediction.match)

table(GSE.query$predicted.id)

Idents(GSE.query) <- "predicted.id"
#GSE.query[["RNA"]] <- JoinLayers(GSE.query[["RNA"]])

saveRDS(GSE.query,"~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE200487_fap_steven_label_transfer.rds")

GSE.query <- readRDS("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE200487_fap_steven_label_transfer.rds")


Idents(GSE.query) <- "adjusted_celltype"
cd55 <- WhichCells(GSE.query, idents = c("CD55"))
dlk1 <- WhichCells(GSE.query, idents = c("DLK1"))
ccn3 <- WhichCells(GSE.query, idents = c("CCN3"))
cd74 <- WhichCells(GSE.query, idents = c("CD74"))
tnmd <- WhichCells(GSE.query, idents = c("TNMD"))
gli1 <- WhichCells(GSE.query, idents = c("GLI1"))
piezo2 <- WhichCells(GSE.query, idents = c("PIEZO2"))
atf3 <- WhichCells(GSE.query, idents = c("ATF3"))
lama2 <- WhichCells(GSE.query, idents = c("LAMA2"))
tnnt1 <- WhichCells(GSE.query, idents = c("TNNT1"))
myl9 <- WhichCells(GSE.query, idents = c("MYL9"))



pdf("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Plot/GSE200487/1_Steven_UMAP.pdf", width = 11, height = 7 )
DimPlot(GSE.query, label=T, group.by="adjusted_celltype", cells.highlight = list(dlk1, cd55, gli1,cd74,lama2,atf3,ccn3, tnmd,piezo2,myl9,tnnt1), cols.highlight = c("#f57e20", "#6a3e97", "#a6cfe3", "#f59999","#1a79b6","#b3d88a","#fabe6e","#e32026","#e9e71f", "#c9b1d3","#32a048"), cols= "grey", pt.size = 0.5, order = T, sizes.highlight = 0.5, label.size = 0)
dev.off()

Idents(GSE.query) <- "dataset"
GSE200487 <- subset(GSE.query, idents = "GSE200487")
Idents(GSE200487) <- "predicted.id"
cd55 <- WhichCells(GSE200487, idents = c("CD55"))
dlk1 <- WhichCells(GSE200487, idents = c("DLK1"))
ccn3 <- WhichCells(GSE200487, idents = c("CCN3"))
cd74 <- WhichCells(GSE200487, idents = c("CD74"))
tnmd <- WhichCells(GSE200487, idents = c("TNMD"))
gli1 <- WhichCells(GSE200487, idents = c("GLI1"))
piezo2 <- WhichCells(GSE200487, idents = c("PIEZO2"))
atf3 <- WhichCells(GSE200487, idents = c("ATF3"))
lama2 <- WhichCells(GSE200487, idents = c("LAMA2"))
tnnt1 <- WhichCells(GSE200487, idents = c("TNNT1"))
myl9 <- WhichCells(GSE200487, idents = c("MYL9"))

pdf("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Plot/GSE200487/1_200487_UMAP.pdf", width = 11, height = 7 )
DimPlot(GSE.query, label=T, group.by="adjusted_celltype", cells.highlight = list(cd55, gli1,cd74,lama2,atf3, tnmd,myl9,tnnt1), cols.highlight = c("#e9e71f", "#c9b1d3", "#a6cfe3", "#f59999","#1a79b6","#b3d88a","#fabe6e","#e32026"), cols= "grey", pt.size = 0.5, order = T, sizes.highlight = 0.5, label.size = 5)
dev.off()


DimPlot(GSE.query, split.by = "dataset", group.by = "predicted.id")

Idents(GSE.query) <- "dataset"
GSE <- subset(GSE.query, idents = "GSE200487")
steven <- subset(GSE.query, idents = "steven")

GSE <- NormalizeData(GSE, normalization.method = "LogNormalize", scale.factor = 10000)
GSE <- FindVariableFeatures(GSE, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE)
GSE <- ScaleData(GSE, features = all.genes)
GSE <- RunPCA(GSE, features = VariableFeatures(object = GSE), npcs = 30)

for (j in 10:12){
  #cluster!
  GSE <- FindNeighbors(GSE, dims = 1:30)
  GSE <- FindClusters(GSE, resolution = (j/10))
  head(Idents(GSE), 5)
  
  #non-linear reduction to umap
  GSE <- RunUMAP(GSE, dims = 1:30)
  DimPlot(GSE, reduction = "umap", label = TRUE)
}



Idents(GSE) <- "predicted.id"
DimPlot(GSE, split.by = c("cluster"))



pdf("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Plot/GSE167186/gse167186_split_feature_plot_2.pdf", width = 10, height = 20)
FeaturePlot(GSE.query, features = c("CD55","DLK1","CD74","DPP4","MME","GLI1"), cols = c("#b1e0e6","#8b0000"),order = T, split.by = "dataset")
dev.off()

FeaturePlot_scCustom(GSE.query, features = c("CD55"), split.by = "dataset", )


test <- readRDS("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Data/GSE200487/FAPS30_all_E1.rds")

# GSE ---------------------------------------------------------------------
#gse on steven
Idents(GSE.query) <- "cluster"
cxcl14 <- WhichCells(GSE.query, idents = c("Cxcl14_FAPs"))
dpp4 <- WhichCells(GSE.query, idents = c("Dpp4_FAPs"))
wisp1 <- WhichCells(GSE.query, idents = c("Wisp1_FAPs"))
dlk1 <- WhichCells(GSE.query, idents = c("Dlk1_FAPs"))
osr1 <- WhichCells(GSE.query, idents = c("Osr1_FAPs"))
activated <- WhichCells(GSE.query, idents = c("ActivatedFAPs"))

pdf("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Plot/GSE138826/2_138826_UMAP.pdf", width = 11, height = 7 )
DimPlot(GSE.query, label=T, group.by="cluster", cells.highlight = list(dlk1, dpp4, cxcl14,wisp1,activated,osr1), cols.highlight = c("#fabe6e", "#b3d88a", "#c9b1d3", "#a6cfe3","#e32026", "#32a048"), cols= "grey", pt.size = 0.5, order = T, sizes.highlight = 0.5,label.size = 0)
dev.off()


Idents(GSE.query) <- "dataset"
steven <- subset(GSE.query, idents = "steven")
Idents(steven) <- "predicted.id"

cxcl14 <- WhichCells(steven, idents = c("Cxcl14_FAPs"))
dpp4 <- WhichCells(steven, idents = c("Dpp4_FAPs"))
wisp1 <- WhichCells(steven, idents = c("Wisp1_FAPs"))
dlk1 <- WhichCells(steven, idents = c("Dlk1_FAPs"))
osr1 <- WhichCells(steven, idents = c("Osr1_FAPs"))
activated <- WhichCells(steven, idents = c("ActivatedFAPs"))

pdf("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Plot/GSE138826/2_Steven_UMAP.pdf", width = 11, height = 7 )
DimPlot(GSE.query, label=T, group.by="cluster", cells.highlight = list(dlk1, dpp4, cxcl14,wisp1,activated,osr1), cols.highlight = c("#fabe6e", "#b3d88a", "#c9b1d3", "#a6cfe3","#e32026", "#32a048"), cols= "grey", pt.size = 0.5, order = T, sizes.highlight = 0.5,label.size = 0)
dev.off()


