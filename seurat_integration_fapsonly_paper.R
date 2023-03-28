suppressPackageStartupMessages({
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(magrittr)
library(Matrix)
})


integrated<-readRDS(file= "../only_faps_integrated.rds")

DefaultAssay(integrated) <- "RNA"


#normalize
integrated <- NormalizeData(integrated, normalization.method = "LogNormalize", scale.factor = 10000)

#identify variable features
integrated <- FindVariableFeatures(integrated, assay = "RNA")

DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs=40)
integrated <- RunUMAP(integrated, reduction="pca",dims=1:30)

seur<-RunTSNE(seur, )


#easier loop res 3->10
for (i in 1:10){
  #cluster!
  integrated <- FindNeighbors(integrated, dims = 1:30)
  integrated <- FindClusters(integrated, resolution = (i/10))
  head(Idents(integrated), 5)
  
  #non-linear reduction to umap
  integrated <- RunUMAP(integrated, dims = 1:30)
}


saveRDS(integrated, file= "../all_faps.rds" )

#find diff expressed genes - all clusters

integrated.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)


#save markers in .csv file
write.csv(integrated.markers, file=paste (".."))


####subset
integrated <- AddMetaData(object = integrated, metadata = integrated@active.ident, col.name = "muscle")
Idents(integrated)<-"orig.ident"
integrated <- RenameIdents(object=integrated, 'seur'="cuff_42")


new.cluster.ids <- c("cuff","delt","cuff",
                     "delt","delt","cuff",
                     "cuff","delt","cuff",
                     "delt","ham","ham",
                     "delt","cuff")
names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)
#

####label partial and full tears
Idents(integrated)<-"orig.ident"
integrated <- AddMetaData(object = integrated, metadata = integrated@active.ident, col.name = "tear")
Idents(integrated)<-"orig.ident"
integrated <- RenameIdents(object=integrated, 'seur'="cuff_42")


new.cluster.ids <- c("partial_cuff","full_delt","full_cuff",
                     "partial_delt","full_delt","full_cuff",
                     "partial_cuff","partial_delt","full_cuff",
                     "full_delt","ham","ham",
                     "partial_delt","partial_cuff")
names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)
integrated$tear<-integrated@active.ident

##rep analyses
prop<-table(Idents(integrated),integrated$orig.ident)
prop<-prop.table(table(Idents(integrated),integrated$orig.ident), margin=2)
prop<-prop.table(table(Idents(integrated),integrated$integrated_snn_res.0.3), margin=2)
write.csv(prop, file=paste ("../prop_table_res3.csv"))


##cuff vs delt by tear type
Idents(integrated)<-"tear"
integrated.markers <- FindMarkers(object=integrated, ident.1="partial_delt",ident.2="partial_cuff", only.pos = F, min.pct = 0.1, logfc.threshold = 0.01)

write.csv(integrated.markers, file=paste ("partialtear_vs_delt.csv"))

integrated.markers <- FindMarkers(object=integrated, ident.1="full_delt",ident.2="full_cuff", only.pos = F, min.pct = 0.1, logfc.threshold = 0.01)

write.csv(integrated.markers, file=paste ("fulltear_vs_delt.csv"))

integrated.markers <- FindMarkers(object=integrated, ident.1="partial_cuff",ident.2="full_cuff", only.pos = F, min.pct = 0.1, logfc.threshold = 0.25)

write.csv(integrated.markers, file=paste ("partialtear_vs_full.csv"))

###now do cuff vs delt for each cluster

#loop though..
Idents(integrated)<-"integrated_snn_res.0.3"

pops<-c("0","1","2","3","4","5","6","7","8","9","10")

for (i in 1:11){
  
  integrated.markers <- FindMarkers(object=integrated, 
                                    subset.ident = pops[i],
                                    group.by="tear",
                                    ident.1="full_delt",
                                    ident.2="full_cuff",
                                    only.pos = F, min.pct = 0.1, logfc.threshold = 0.25)
  write.csv(integrated.markers, file=paste0("../fulldelt_vs_cuff_cluster",pops[i],".csv"))
  
}


#make images for figure 1
library(RColorBrewer)
DefaultAssay(integrated)<- "RNA"
Idents(integrated)<-"integrated_snn_res.0.3"
DimPlot(integrated, reduction = "umap", cols = brewer.pal(11, "Paired"), label = F, pt.size = 1, raster = T, raster.dpi = c(600,600))+NoLegend()+theme_void()

gene<-c("CCN3","LAMA2","CD74","DLK1","ATF3","CD55","GLI1","PIEZO2","MYL9","TNMD","TNNT1")

FeaturePlot(integrated, features = gene, pt.size = 2, order = T, cols = c("powderblue","red4"), raster = TRUE, min.cutoff = "q10",ncol=4,raster.dpi = c(600,600))&theme_void()&NoLegend()

DotPlot(integrated, group.by = "cells",features = gene,cols = c("powderblue","red4"),dot.scale = 8, col.min=0, dot.min=.10) + 
  RotatedAxis()

##make images for sup fig 1

Idents(integrated)<-"integrated_snn_res.0.3"
integrated$order<-factor(integrated$orig.ident, 
                  levels=c("delt_42","seur",
                           "delt_43m_1.25.22","cuff_43m_1.25.22",
                           "delt_48f_2.10.22","cuff_48f_2.10.22",
                           "delt_72_1.11.22","cuff_72_1.11.22",
                           "delt_56m_1.19.22","cuff_56m_1.19.22",
                           "delt_79f_1.25.22","cuff_79f_1.25.22",
                           "twentyfour","thirtytwo"))

DimPlot(integrated, reduction = "umap", cols = brewer.pal(11, "Paired"), label = F, pt.size = 5, raster = T, raster.dpi = c(600,600),split.by="order",ncol=2)&NoLegend()&theme_void()

Idents(integrated)<-"integrated_snn_res.0.3"
integrated.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.25)
integrated <- ScaleData(integrated)


integrated.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

top5 = integrated.markers %>% group_by(cluster) %>% slice_max(order_by=avg_log2FC, n=5)

DoHeatmap(integrated, features =unique(top5$gene),assay = "RNA",raster = T,group.by="integrated_snn_res.0.3",group.colors = brewer.pal(11, "Paired")) 


#make images for figure 2 and sup figure 2

library(RColorBrewer)
DefaultAssay(integrated)<- "RNA"
Idents(integrated)<-"integrated_snn_res.0.3"
p<-c("#6BAED6","#08519C","#2171B5","#A1D99B","#FDAE6B","#CB181D",
     "#F16913","#4292C6","#9ECAE1","#54278F","#A50F15","#807DBA","#EF3B2C")
DimPlot(integrated, reduction = "umap", cols = "Paired", label = F, pt.size = 1, raster = T, raster.dpi = c(600,600),split.by="muscle")+NoLegend()+theme_void()


pdf(file=".pdf", 4.8, 4.9)
FeaturePlot(integrated, features = "OSR1", pt.size = 1, order = T, cols = c("powderblue","red4"), raster = TRUE, max.cutoff="q90",raster.dpi = c(600,600)) +theme_void()+NoLegend()
dev.off()


#sup figure 2

gene<-c("CCN3","PENK","PI16","THBS4")
gene<-c("LAMA2","SCN7A","ITIH5","LINC00324")
gene<-c("ACTA2","MUSTN1","ATF3","EGR1")
gene<-c("MYL9","TAGLN2","CD248","RHOC")
gene<-c("TNMD","POSTN","IGFBP2","ITGA10")
gene<-c("TNNT1","ACTA1","TNNC2","DES")


pdf(file="...pdf", 4.8, 4.9)
FeaturePlot(integrated, features = gene, pt.size = 1, order = T, cols = c("powderblue","red4"), raster = TRUE, max.cutoff="q90",raster.dpi = c(600,600))&theme_void()&NoLegend()
dev.off()

#make images for Figure 5

#PRESTUFF
Idents(integrated)<-"tear"
integrated.markers <- FindMarkers(object=integrated, ident.1="partial_delt",ident.2="partial_cuff", only.pos = F, min.pct = 0.01, logfc.threshold = 0)

write.csv(integrated.markers, file=paste ("../partialtear_vs_delt.csv"))

integrated.markers <- FindMarkers(object=integrated, ident.1="full_delt",ident.2="full_cuff", only.pos = F, min.pct = 0.01, logfc.threshold = 0)

write.csv(integrated.markers, file=paste ("../fulltear_vs_delt.csv"))

#


partial.markers<- read.csv(file=paste ("partialtear_vs_delt.csv"), row.names = "X")
full.markers<- read.csv(file=paste ("fulltear_vs_delt.csv"), row.names = "X")

full.markers$gene<-rownames(full.markers)
partial.markers$gene<-rownames(partial.markers)

int<-intersect(rownames(partial.markers),rownames(full.markers))

int_sub_full<-full.markers[rownames(full.markers) %in% int,]
int_sub_partial<-partial.markers[rownames(partial.markers) %in% int,]

int_sub_full<-data.frame(gene=rownames(int_sub_full),log2FC=int_sub_full$avg_log2FC)
int_sub_partial<-data.frame(gene=rownames(int_sub_partial),log2FC=int_sub_partial$avg_log2FC)


all.merge.v2<-merge(int_sub_full,int_sub_partial,by="gene",all.x=TRUE,all.y=TRUE)
colnames(all.merge.v2)[2]<-"Full.Tear"
colnames(all.merge.v2)[3]<-"Partial.Tear"


#remove genes if Log2FC not >.322 in at least one condition
all.merge.v2=subset(all.merge.v2[all.merge.v2$Full.Tear<(-.1) | 
                                         all.merge.v2$Full.Tear>.1 |
                                   all.merge.v2$Partial.Tear<(-.1) | 
                                      all.merge.v2$Partial.Tear>.1                                     ,])
#


#label points
genel<-c("THY1","PLA2G2A","COL14A1","CD55","CCN3","CLU","TPPP3","FN1",
         "FGF18","FBN1","CD34","MEG3","MEG8","COL15A1","ITGBL1",
         "FOS","FOSB","CCN1","ACTA1","CCN3","ADH1B","MYC",
         "CDKN1A","CCN2","CCL2","PTGDS","ITGA8","CXCL14",
         "COL4A1","EGFR","LAMA2","PI16",
         "PDGFRA","PDE5A","SOX9","OSR2","OSR1",
         "COL12A1","GDF10","DPP4","DLK1","MME")  
genel<-c("THY1","PLA2G2A","COL14A1","CD55","CCN3","CLU","TPPP3","FN1",
         "FGF18","CD34","COL15A1",
         "FOS","FOSB","CCN1","ACTA1","CCN3","ADH1B","MYC",
         "CDKN1A","CCN2","CCL2","PTGDS","ITGA8","CXCL14",
         "LAMA2","PI16",
         "PDGFRA","SOX9","OSR2","OSR1",
         "GDF10","DPP4","DLK1","MME") 
cdat <- all.merge.v2 %>%
  mutate(plotname = as.character(gene))
cdat <- cdat %>%
  mutate(plotname = ifelse(plotname %in% genel, plotname, ""))
cdat <- cdat %>%
  mutate(colorit = ifelse(plotname %in% genel, 1, 0))

highlight_df <- cdat %>% 
  filter(colorit>0)

library(ggrastr)
library(ggrepel)

p<-ggplot(cdat,aes(x=Full.Tear,y= Partial.Tear))+
  rasterise(geom_point(size=1,color="grey"))+
  ggtitle("Partial and Full Tears") + 
  xlab("Full Tear (Log2FC)") + ylab("Partial Tear (Log2FC)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0))+
  geom_hline(yintercept = 0.32, linetype="dashed",color="blue")+
  geom_hline(yintercept = -0.32, linetype="dashed",color="blue")+
  geom_vline(xintercept = 0.32, linetype="dashed",color="blue")+
  geom_vline(xintercept = -0.32, linetype="dashed",color="blue")+
  geom_label_repel(aes(label=plotname),
                  max.overlaps = Inf,
                  box.padding=0.5,
                  min.segment.length = 0,
                  size=3,
                  segment.size=0.25,
                  segment.colour="black",
                  segment.linetype=2,
                  force=50)+
  geom_point(data=highlight_df,
             aes(x=Full.Tear,y= Partial.Tear),
             color='black',
             size=2)+
  ylim(3,-3)+
  xlim(3,-3)

#volcano
library(EnhancedVolcano)

full.markers<- read.csv(file=paste ("fulltear_vs_delt.csv"))


gene<-c("THY1","DLK1","CD55","CCN3","COL4A1","PRG4","SEMA3C",
        "CXCL14","FGF18","DPP4","COL14A1",
        "ACTA1","PLA2G2A","ADH1B","MME","SOX9",
        "COL15A1","MEG3","CLEC14A","LAMA2",
        "PDGFRA","FOXO3","PTGDS")

pdf(file="..", 7, 5)
EnhancedVolcano(full.markers,
                lab = full.markers$X,
                selectLab = gene,
                boxedLabels = F,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Delt vs Cuff',
                pCutoff = 10e-12,
                FCcutoff = 0.322,
                pointSize = c(ifelse(full.markers$X %in% gene,3,1)),
                labSize = 5.0, 
                colAlpha = 4/5,
                legendPosition = 'none',
                legendLabSize = 24,
                legendIconSize = 4.0,
                gridlines.major = F,
                gridlines.minor = F,
                drawConnectors = T,
                widthConnectors = 0.5,
                raster = T)+
  ggplot2::coord_cartesian(xlim=c(-3, 3))+
  scale_x_reverse()+
  scale_fill_continuous(guide=FALSE)+
  coord_flip()
dev.off()

# make images for Figure 6 and sup 6


Idents(integrated)<-"tear"
integrated=subset(integrated, idents=c("full_delt","full_cuff"))
Idents(integrated)<-"integrated_snn_res.0.3"


gene<-c("")

for (i in 1:length(gene)){
pdf(file=paste0("fulltear_vs_delt_ident3_raster_5x5_600dpi",gene[i],".pdf"), 5,3)
print(VlnPlot(integrated,idents = "3", split.by = "tear",features = gene[i],raster = T,pt.size = 0,flip=T,cols = c("powderblue","powderblue","red4","red4","red4"), fill.by = "ident")+ theme_void()+NoLegend())
dev.off()
}


#sup fig 7

gene<-c("CXCL12")

DotPlot(integrated, group.by = "cells",features = gene,cols = c("powderblue","red4"), dot.scale = 8, col.min=0, dot.min=.10) + 
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

Idents(integrated)<-"integrated_snn_res.0.3"
integrated <- RenameIdents(object=integrated, 
                           "0"="CCN3+",
                           "1"="LAMA2+",
                           "2"="CD74+",
                           "3"="DLK1+",
                           "4"="ATF3+",
                           "5"="CD55+",
                           "6"="GLI1+",
                           "7"="PIEZO2+",
                           "8"="MYL9+",
                           "9"="TNMD+",
                           "10"="TNsNT1+")
integrated <- AddMetaData(object = integrated, metadata = integrated@active.ident, col.name = "cells")
Idents(integrated)<-"cells"

integrated <- AddMetaData(object = integrated, metadata = Embeddings(integrated, reduction = "umap")[,1], col.name = "UMAP_1")
integrated <- AddMetaData(object = integrated, metadata = Embeddings(integrated, reduction = "umap")[,2], col.name = "UMAP_2")

write.csv(subset(integrated@meta.data, select = c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","integrated_snn_res.0.3","cells","muscle","sample","tear","UMAP_1","UMAP_2")),
          file = "../fapsonly_metadata_paper.csv")

writeMM(integrated@assays$RNA@counts, file = "../fapsonly_counts_paper.mtx")

