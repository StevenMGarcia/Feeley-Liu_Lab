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
    #library(nichenetr)
  })

data_path <- "~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Data/"
mouse_data <- readRDS(paste0(data_path,"GSE138826/GSE138826_regen_data.rds"))

mouse <- CreateSeuratObject(counts = mouse_data$RNA@data, meta.data = mouse_data$RNA@varMetadata)

obj <- mouse

Idents(obj) <- "metacluster"
  # QC ----------------------------------------------------------------------
#  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
#  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
#  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
#  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#  plot1 + plot2
  
 # obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
 # feature_RNA <- 0.83*max(obj@meta.data$nFeature_RNA)
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
  
  Idents(obj) <- "metacluster"
  DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1) #preliminary clustering
  
# Save merged object ------------------------------------------------------
saveRDS(obj, file = "~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE138826.rds")
  
GSE138826 <- obj
rm(obj)

# Load Cells ------------------------------------------------------------------
GSE138826 <- readRDS("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE138826.rds")
GSE138826@meta.data$dataset <- "GSE138826"
Idents(GSE138826) <- "metacluster"


GSE_fap <- subset(GSE138826, idents = "FAPs")
# Mouse to Human  ----------------------------------------------------------
exp_mtx <- as.matrix(GSE_fap@assays$RNA$data)

con_df <- data.frame(hum_orig = rownames(exp_mtx),
                     mouse = convert_mouse_to_human_symbols(rownames(exp_mtx)),
                     stringsAsFactors = F)
head(con_df,10)
con_df <- con_df[!is.na(con_df$mouse),,F]
head(con_df,10)

exp_mtx <- exp_mtx[con_df$hum_orig,]
rownames(exp_mtx) <- make.names(con_df$mouse, unique = T)


gse_hum_fap <- CreateSeuratObject(counts = exp_mtx, meta.data = GSE_fap@meta.data)

gse_hum_fap <- SCTransform(gse_hum_fap)
# These are now standard steps in the Seurat workflow for visualization and clustering
gse_hum_fap <- RunPCA(gse_hum_fap, verbose = FALSE)
gse_hum_fap <- RunUMAP(gse_hum_fap, dims = 1:30, verbose = FALSE)

gse_hum_fap <- FindNeighbors(gse_hum_fap, dims = 1:30, verbose = FALSE)
gse_hum_fap <- FindClusters(gse_hum_fap, verbose = FALSE)
DimPlot(gse_hum_fap, label = TRUE) #preliminary clustering

saveRDS(gse_hum_fap, file = "~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE138826_fap_converted_to_human.rds")

gse_hum_fap <- readRDS("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE138826_fap_converted_to_human.rds")
# Load Steven -------------------------------------------------------------
steven_fap <- readRDS("~/Desktop/UCSF/Seurat_Objects/SMG_FAP.rds")
#steven_fap@meta.data$dataset <- "steven"
#DimPlot(steven_fap)
#steven_fap@meta.data$adjusted_celltype <- steven_fap@meta.data$celltype

#levels(steven_fap@meta.data$adjusted_celltype) <- c("CCN3","LAMA2","CD74","DLK1","ATF3","CD55","GLI1","PIEZO2","MYL9","TNMD","TNNT1")
Idents(steven_fap) <- "adjusted_celltype"
DimPlot(steven_fap)
#saveRDS(steven_fap, "~/Desktop/UCSF/Seurat_Objects/SMG_FAP.rds")

# Merge Cells -------------------------------------------------------------
merged <- merge(steven_fap,y = gse_hum_fap, add.cell.ids = c("steven","GSE138826"))

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



saveRDS(merged, "/~Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE138826_fap_converted_to_human_steven.rds")
merged <- readRDS("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE138826_fap_converted_to_human_steven.rds")


# Label Transfer ----------------------------------------------------------
#Idents(merged) <- "cluster"
#Idents(merged) <- "dataset"

steven.ref <- steven_fap
GSE.query <- gse_hum_fap
#rm(merged)

DefaultAssay(GSE.query) <- "RNA"
DefaultAssay(steven.ref) <- "RNA"


GSE.query <- NormalizeData(GSE.query)


cell.anchors <- FindTransferAnchors(reference = steven.ref, 
                                    query = GSE.query, 
                                    dims = 1:40, 
                                    reference.reduction = "pca")

predictions <- TransferData(anchorset = cell.anchors, 
                            refdata = steven.ref$adjusted_celltype, 
                            dims = 1:25)
GSE.query <- AddMetaData(GSE.query, metadata = predictions)
Idents(GSE.query) <- "predicted.id"

#GSE.query[["RNA"]] <- JoinLayers(GSE.query[["RNA"]])

DimPlot_scCustom(GSE.query, label = F,figure_plot = T, colors_use = 
c("#a6cfe3","#e32026","#c9b1d3","#f59999","#f57e20","#32a048","#1a79b6","#fabe6e","#b3d88a","#e9e71f","#6a3e97"))

#DimPlot_scCustom(GSE.query, split.by = "dataset", label = F,figure_plot = T, colors_use = c("#e32026", "#a6cfe3", "#6a3e97", "#32a048","#1a79b6","#f57e20"))

saveRDS(GSE.query,"Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE138826_fap_steven_label_transfer.rds" )

GSE.query = readRDS("Desktop/UCSF/Steven/HFAP_Integration_other_Data/Objects/GSE138826_fap_steven_label_transfer.rds")

#STeven on other
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



pdf("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Plot/GSE138826/1_Steven_UMAP.pdf", width = 11, height = 7 )
DimPlot(GSE.query, label=T, group.by="adjusted_celltype", cells.highlight = list(dlk1, cd55, gli1,cd74,lama2,atf3,ccn3, tnmd,piezo2,myl9,tnnt1), cols.highlight = c("#f57e20", "#6a3e97", "#a6cfe3", "#f59999","#1a79b6","#b3d88a","#fabe6e","#e32026","#e9e71f", "#c9b1d3","#32a048"), cols= "grey", pt.size = 0.5, order = T, sizes.highlight = 0.5, label.size = 0)
dev.off()

Idents(GSE.query) <- "dataset"
GSE138826 <- subset(GSE.query, idents = "GSE138826")
Idents(GSE138826) <- "predicted.id"
cd55 <- WhichCells(GSE138826, idents = c("CD55"))
dlk1 <- WhichCells(GSE138826, idents = c("DLK1"))
ccn3 <- WhichCells(GSE138826, idents = c("CCN3"))
cd74 <- WhichCells(GSE138826, idents = c("CD74"))
tnmd <- WhichCells(GSE138826, idents = c("TNMD"))
gli1 <- WhichCells(GSE138826, idents = c("GLI1"))
piezo2 <- WhichCells(GSE138826, idents = c("PIEZO2"))
atf3 <- WhichCells(GSE138826, idents = c("ATF3"))
lama2 <- WhichCells(GSE138826, idents = c("LAMA2"))
tnnt1 <- WhichCells(GSE138826, idents = c("TNNT1"))
myl9 <- WhichCells(GSE138826, idents = c("MYL9"))

pdf("~/Desktop/UCSF/Steven/HFAP_Integration_other_Data/Plot/GSE138826/1_138826_UMAP.pdf", width = 11, height = 7 )
DimPlot(GSE.query, label=T, group.by="adjusted_celltype", cells.highlight = list(dlk1, cd55, gli1,cd74,lama2,atf3,ccn3, tnmd,piezo2,myl9,tnnt1), cols.highlight = c("#f57e20", "#6a3e97", "#a6cfe3", "#f59999","#1a79b6","#b3d88a","#fabe6e","#e32026","#e9e71f", "#c9b1d3","#32a048"), cols= "grey", pt.size = 0.5, order = T, sizes.highlight = 0.5, label.size = 0)
dev.off()



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








#Integrated UMAP
DimPlot_scCustom(GSE.query, group.by = "dataset", label = F, colors_use = c("firebrick", "cadetblue2"),figure_plot = T)

#Split UMAP with Label Transfer
#DimPlot_scCustom(GSE.query, split.by = "dataset", label = F,figure_plot = T, colors_use = c("#b3d88a", "#a6cfe3", "#6a3e97", "#32a048","#fabe6e","#f57e20","#f59999","#e32026","#1a79b6", "#c9b1d3","#e9e71f"))

Idents(GSE.query) <- "dataset"
GSE <- subset(GSE.query, idents = "GSE138826")
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

FeaturePlot_scCustom(GSE.query, features = c("CD55","DLK1","CD74","DPP4","MME","GLI1"), colors_use = c("cadetblue1","firebrick4"),order = T, num_columns = 3)







