## DF DF DF ##

library(DoubletFinder)
library(Seurat)
library(Matrix)
library(fields)
library(KernSmooth)
library(ROCR)
library(parallel)

#
seur = readRDS("../Adipo.rds")
Idents(seur) <- "seurat_clusters"
DimPlot(seur, reduction = "umap", label = TRUE, pt.size = 1)

## pK Identification (no ground-truth) 

sweep.res.list_seur <- paramSweep_v3(seur, PCs = 1:10, sct = FALSE)
sweep.stats_seur <- summarizeSweep(sweep.res.list_seur, GT = FALSE)
bcmvn_seur <- find.pK(sweep.stats_seur)

## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(seur@meta.data$seurat_clusters)    
nExp_poi <- round(0.075*nrow(seur@meta.data))  ## Assuming 7.5% doublet formation rate - tailor 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur <- doubletFinder_v3(seur, PCs = 1:10, pN = 0.25, pK = 0.12, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

seur <- doubletFinder_v3(seur, PCs = 1:10, pN = 0.25, pK = 0.12, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_279", sct = FALSE)

DimPlot(seur, reduction = "umap", label = TRUE, pt.size = 1, group.by = "pANN_0.25_0.09_279")
DimPlot(seur, reduction = "umap", label = TRUE, pt.size = 1, group.by = "DF.classifications_0.25_0.12_279")


