suppressPackageStartupMessages({
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
})

integrated = readRDS("../Std_Adipo_TGFB_v1.rds")

#sup figure 4
DefaultAssay(integrated)<- "RNA"
Idents(integrated)<-"integrated_snn_res.0.3"


p<-c("#6BAED6","#08519C","#2171B5","#A1D99B","#FDAE6B","#CB181D",
     "#F16913","#4292C6","#9ECAE1","#54278F","#A50F15","#807DBA","#EF3B2C")
DimPlot(integrated, reduction = "umap", cols = "Paired", label = F, pt.size = 1, split.by="orig.ident",raster = T, raster.dpi = c(600,600))+NoLegend()+theme_void()

DimPlot(integrated, reduction = "umap", label = F, pt.size = 1, group.by="orig.ident",raster = T, raster.dpi = c(600,600))+NoLegend()+theme_void()

DimPlot(integrated, reduction = "umap", cols = "Paired", label = F, pt.size = 1, raster = T, group.by="integrated_snn_res.0.3",raster.dpi = c(600,600))+NoLegend()+theme_void()

gene <- c("PDGFRA","CD34","OSR1","OSR2",
          "THY1","GLI1","CD74","DLK1",
          "MME","CD55","DPP4","PIEZO2",
          "SCX","PPARG","PLIN1","ACTA2")


FeaturePlot(integrated, features = gene, pt.size = 1, order = T, cols = c("powderblue","red4"), raster = TRUE, max.cutoff="q90",raster.dpi = c(600,600))&theme_void()&NoLegend()

Idents(integrated)<-"integrated_snn_res.0.3"
integrated$orig.ident<-factor(integrated$orig.ident, 
                         levels=c("std","Adipo","TGFB"))

VlnPlot(integrated, features = gene,raster = F,pt.size = 0, adjust=1, same.y.lims = F, stack = T, fill.by = "ident", flip=T, sort=F,
        cols=c("#619CFF","#F8766D","#00BA38"))


Idents(integrated)<-"orig.ident"
integrated.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)

write.csv(integrated.markers, file=paste ("../Std_Adipo_fibro_markers_res3.csv"))

Idents(integrated)<-"orig.ident"
integrated.markers <- FindAllMarkers(integrated,only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)

write.csv(integrated.markers, file=paste ("../Std_Adipo_fibro_markers_origident.csv"))


#save metadata

integrated <- AddMetaData(object = integrated, metadata = Embeddings(integrated, reduction = "umap")[,1], col.name = "UMAP_1")
integrated <- AddMetaData(object = integrated, metadata = Embeddings(integrated, reduction = "umap")[,2], col.name = "UMAP_2")

write.csv(subset(integrated@meta.data, select = c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","integrated_snn_res.0.7","UMAP_1","UMAP_1")),
          file = "../Std_Adipo_TGFB_metadata_paper.csv")

writeMM(integrated@assays$RNA@counts, file = "../Std_Adipo_TGFB_counts_paper.mtx")
 