suppressPackageStartupMessages({
library(Seurat)
library(dplyr)
library(ggplot2)
})

integrated <- readRDS(file = "../allcells_integrated.rds" )
DefaultAssay(integrated)<- "RNA"

Idents(integrated)<-"integrated_snn_res.0.5"
DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 1, raster = TRUE)


#find diff expressed genes - all clusters
integrated.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.25)

#save markers in .csv file
write.csv(integrated.markers, file=paste ("../allcells_res5.csv"))


#-----------------

DefaultAssay(integrated)<- "RNA"

####subset faps

DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 1)
integrated<-subset(integrated, idents = c("0","1","2","4","7","8","11"))
DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 1)

saveRDS(integrated, file = "../only_faps_integrated.rds")



#make images for figure 1
library(RColorBrewer)
DefaultAssay(integrated)<- "RNA"
Idents(integrated)<-"integrated_snn_res.0.5"
p<-c("deepskyblue","deepskyblue1","deepskyblue2","snow3",
     "deepskyblue3","firebrick3",
     "green4","deepskyblue4","dodgerblue","darkorange2",
     "snow4","dodgerblue3","darkorchid","firebrick","darkorange",
     "darkorchid1","firebrick4","red3")
DimPlot(integrated, reduction = "umap", cols = p, label = F, pt.size = 1, raster = T, raster.dpi = c(600,600))+NoLegend()+theme_void()

#sup fig 1
integrated$order<-factor(integrated$orig.ident, 
                         levels=c("delt_42","seur",
                                  "delt_43m_1.25.22","cuff_43m_1.25.22",
                                  "delt_48f_2.10.22","cuff_48f_2.10.22",
                                  "delt_72_1.11.22","cuff_72_1.11.22",
                                  "delt_56m_1.19.22","cuff_56m_1.19.22",
                                  "delt_79f_1.25.22","cuff_79f_1.25.22",
                                  "twentyfour","thirtytwo"))
DimPlot(integrated, reduction = "umap", cols = p, label = F, pt.size = 5, raster = T, raster.dpi = c(600,600), split.by="order", ncol = 2)+NoLegend()+theme_void()

gene<-c("PDGFRA","CD34","PTPRC","CD68","CD3E","PECAM1","RGS5","ACTA2","PAX7","MYOD1")

FeaturePlot(integrated, features = gene, pt.size = 1, order = T, cols = c("powderblue","red4"), raster = TRUE,ncol=5, raster.dpi = c(600,600))&theme_void()&NoLegend()

#sup figure 7
Idents(integrated)<-"integrated_snn_res.0.5"
gene<-c("CXCR4","PTPRC","PAX7","PDGFRA","CD68","CD3E","CXCL12","DHH","PECAM1")

DotPlot(integrated, features = gene,cols = c("powderblue","red4")) + 
  RotatedAxis()


#save metadata file

Idents(integrated)<-"orig.ident"
integrated <- RenameIdents(object=integrated, 
                           "cuff_43m_1.25.22"="partial_cuff_43M",
                           "cuff_48f_2.10.22"="partial_cuff_48F",
                           "cuff_56m_1.19.22"="full_cuff_56M",
                           "cuff_72_1.11.22"="full_cuff_72F",
                           "cuff_79f_1.25.22"="full_cuff_79F",
                           "delt_42"="partial_delt_42F",
                           "delt_43m_1.25.22"="partial_delt_43M",
                           "delt_48f_2.10.22"="partial_delt_48F",
                           "delt_56m_1.19.22"="full_delt_56M",
                           "delt_72_1.11.22"="full_delt_72F",
                           "delt_79f_1.25.22"="full_delt_79F",
                           "seur"="partial_cuff_42F",
                           "thirtytwo"="ham32",
                           "twentyfour"="ham24")

integrated <- AddMetaData(object = integrated, metadata = integrated@active.ident, col.name = "sample")
Idents(integrated)<-"sample"

integrated <- RenameIdents(object=integrated, 
                           "partial_cuff_43M"="partial_cuff",
                           "partial_cuff_48F"="partial_cuff",
                           "full_cuff_56M"="full_cuff",
                           "full_cuff_72F"="full_cuff",
                           "full_cuff_79F"="full_cuff",
                           "partial_delt_42F"="partial_delt",
                           "partial_delt_43M"="partial_delt",
                           "partial_delt_48F"="partial_delt",
                           "full_delt_56M"="full_delt",
                           "full_delt_72F"="full_delt",
                           "full_delt_79F"="full_delt",
                           "partial_cuff_42F"="partial_cuff",
                           "ham32"="ham",
                           "ham24"="ham")

integrated <- AddMetaData(object = integrated, metadata = integrated@active.ident, col.name = "tear")
Idents(integrated)<-"tear"


Idents(integrated)<-"integrated_snn_res.0.5"
integrated <- RenameIdents(object=integrated, 
                           "0"="FAP",
                           "1"="FAP",
                           "2"="FAP",
                           "3"="Fibroblast",
                           "4"="FAP",
                           "5"="Endothelial",
                           "6"="Myogenic",
                           "7"="FAP",
                           "8"="FAP",
                           "9"="Pericyte",
                           "10"="Fibroblast",
                           "11"="FAP",
                           "12"="Macrophage",
                           "13"="Endothelial",
                           "14"="Pericyte",
                           "15"="T-cell",
                           "16"="Endothelial",
                           "17"="Endothelial")

integrated <- AddMetaData(object = integrated, metadata = integrated@active.ident, col.name = "cells")
Idents(integrated)<-"cells"



Idents(integrated)<-"orig.ident"

integrated <- AddMetaData(object = integrated, metadata = Embeddings(integrated, reduction = "umap")[,1], col.name = "UMAP_1")
integrated <- AddMetaData(object = integrated, metadata = Embeddings(integrated, reduction = "umap")[,2], col.name = "UMAP_2")

write.csv(subset(integrated@meta.data, select = c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","integrated_snn_res.0.5","cells","muscle","sample","tear","UMAP_1","UMAP_2")),
          file = "../allcells_metadata_paper.csv")

writeMM(integrated@assays$RNA@counts, file = "../allcells_counts_paper.mtx")

