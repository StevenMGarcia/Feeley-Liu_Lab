remove.packages("CellChat")
devtools::install_github("sqjin/CellChat")
# Libraries ---------------------------------------------------------------
library(Seurat)
library(CellChat)
library(patchwork)
library(plyr)
library(NMF)
library(ggalluvial)
library(tidyr)
library(stringr)
# Load Object -------------------------------------------------------------
integrated <- readRDS("~/Desktop/UCSF/Steven/CellChat/hum_muscle_hetero_faps.rds") 

levels(integrated@meta.data$reanalyzed) <-c("FAP 0","FAP 1","CD74","DLK1","FAP 4","CD55","GLI-1","PIEZO2","FAP 8","SCX/ TNMD","FAP 10", "Fibrogenic","Pericytes","Endothelial", "Satellite", "Macrophages", "T-Cells") #new markers
integrated <- SetIdent(integrated, value = "reanalyzed")

integrated <- subset(integrated, idents = c("CD74","DLK1","CD55","GLI-1","PIEZO2","SCX/ TNMD" ,"Fibrogenic","Pericytes","Endothelial", "Satellite", "Fibrogenic", "Macrophages", "T-Cells"))

Idents(integrated)<-"orig.ident"
integrated <- AddMetaData(object = integrated, metadata = integrated@active.ident, col.name = "tear")
Idents(integrated)<-"orig.ident"
integrated <- RenameIdents(object=integrated, 
                           'seur'="cuff_42")
new.cluster.ids <- c("partial_cuff","full_delt","full_cuff",
                     "partial_delt","full_delt","full_cuff",
                     "partial_cuff","partial_delt","full_cuff",
                     "full_delt","ham","ham",
                     "partial_delt","partial_cuff")
names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)
integrated$tear<-integrated@active.ident

integrated <- SetIdent(integrated, value = "reanalyzed")

# tear type ---------------------------------------------------------------
partial_cuff <- subset(integrated, subset = tear == "partial_cuff")
full_cuff <- subset(integrated, subset = tear == "full_cuff")

partial_delt <- subset(integrated, subset = tear == "partial_delt")
full_delt <- subset(integrated, subset = tear == "full_delt")

# pc ----------------------------------------------------------------------
data.input_pc <- GetAssayData(partial_cuff, assay = "RNA", slot = "data") # normalized data matrix
labels_pc <- Idents(partial_cuff)
meta_pc <- data.frame(group = labels_pc, row.names = names(labels_pc)) # create a dataframe of the cell labels
names(meta_pc)[1] <- "labels"

cellchat_pc <- createCellChat(object = data.input_pc, meta = meta_pc, group.by = "labels")

cellchat_pc <- addMeta(cellchat_pc, meta = meta_pc, meta.name = "labels")
cellchat_pc <- setIdent(cellchat_pc, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_pc@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_pc@idents))

CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat_pc@DB <- CellChatDB.use

cellchat_pc <- subsetData(cellchat_pc) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4)

cellchat_pc <- identifyOverExpressedGenes(cellchat_pc)
cellchat_pc <- identifyOverExpressedInteractions(cellchat_pc)  #LR

cellchat_pc <- projectData(cellchat_pc, PPI.human)

cellchat_pc <- computeCommunProb(cellchat_pc) #net
cellchat_pc <- filterCommunication(cellchat_pc, min.cells = 10)

cellchat_pc <- computeCommunProbPathway(cellchat_pc) #netP
cellchat_pc <- aggregateNet(cellchat_pc)

mat_pc <- cellchat_pc@net$weight

#subset fap only communication
#cellchat_pc@netP$prob[7:12,7:12,]<- 0
#cellchat_pc@netP$prob[7:12,7:12,]<- 0 
#cellchat_pc@net$pval[7:12,7:12,]<- 0
#cellchat_pc@net$count[7:12,7:12]<- 0 
#cellchat_pc@net$weight[7:12,7:12]<- 0



# pd ----------------------------------------------------------------------
data.input_pd <- GetAssayData(partial_delt, assay = "RNA", slot = "data") # normalized data matrix
labels_pd <- Idents(partial_delt)
meta_pd <- data.frame(group = labels_pd, row.names = names(labels_pd)) # create a dataframe of the cell labels
names(meta_pd)[1] <- "labels"

cellchat_pd <- createCellChat(object = data.input_pd, meta = meta_pd, group.by = "labels")

cellchat_pd <- addMeta(cellchat_pd, meta = meta_pd, meta.name = "labels")
cellchat_pd <- setIdent(cellchat_pd, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_pd@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_pd@idents))

CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat_pd@DB <- CellChatDB.use

cellchat_pd <- subsetData(cellchat_pd) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4)

cellchat_pd <- identifyOverExpressedGenes(cellchat_pd)
cellchat_pd <- identifyOverExpressedInteractions(cellchat_pd)  #LR

cellchat_pd <- computeCommunProb(cellchat_pd) #net
cellchat_pd <- filterCommunication(cellchat_pd, min.cells = 10)


cellchat_pd <- computeCommunProbPathway(cellchat_pd) #netP
cellchat_pd <- aggregateNet(cellchat_pd)

mat_pd <- cellchat_pd@net$weight
#subset fap only communication
#cellchat_pd@netP$prob[7:12,7:12,]<- 0
#cellchat_pd@netP$prob[7:12,7:12,]<- 0 
#cellchat_pd@net$pval[7:12,7:12,]<- 0
#cellchat_pd@net$count[7:12,7:12]<- 0 
#cellchat_pd@net$weight[7:12,7:12]<- 0

# fc ----------------------------------------------------------------------
data.input_fc <- GetAssayData(full_cuff, assay = "RNA", slot = "data") # normalized data matrix
labels_fc <- Idents(full_cuff)
meta_fc <- data.frame(group = labels_fc, row.names = names(labels_fc)) # create a dataframe of the cell labels
names(meta_fc)[1] <- "labels"

cellchat_fc <- createCellChat(object = data.input_fc, meta = meta_fc, group.by = "labels")

cellchat_fc <- addMeta(cellchat_fc, meta = meta_fc, meta.name = "labels")
cellchat_fc <- setIdent(cellchat_fc, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_fc@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_fc@idents))

CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat_fc@DB <- CellChatDB.use

cellchat_fc <- subsetData(cellchat_fc) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4)

cellchat_fc <- identifyOverExpressedGenes(cellchat_fc)
cellchat_fc <- identifyOverExpressedInteractions(cellchat_fc)  #LR

cellchat_fc <- computeCommunProb(cellchat_fc) #net
cellchat_fc <- filterCommunication(cellchat_fc, min.cells = 10)


cellchat_fc <- computeCommunProbPathway(cellchat_fc) #netP
cellchat_fc <- aggregateNet(cellchat_fc)

mat_fc <- cellchat_fc@net$weight
#subset fap only communication
#cellchat_fc@netP$prob[7:12,7:12,]<- 0
#cellchat_fc@net$prob[7:12,7:12,] <- 0 
#cellchat_fc@net$pval[7:12,7:12,] <- 0
#cellchat_fc@net$count[7:12,7:12,] <- 0 
#cellchat_fc@net$weight[7:12,7:12,] <- 0

# fd ----------------------------------------------------------------------
data.input_fd <- GetAssayData(full_delt, assay = "RNA", slot = "data") # normalized data matrix
labels_fd <- Idents(full_delt)
meta_fd <- data.frame(group = labels_fd, row.names = names(labels_fd)) # create a dataframe of the cell labels
names(meta_fd)[1] <- "labels"

cellchat_fd <- createCellChat(object = data.input_fd, meta = meta_fd, group.by = "labels")

cellchat_fd <- addMeta(cellchat_fd, meta = meta_fd, meta.name = "labels")
cellchat_fd <- setIdent(cellchat_fd, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_fd@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_fd@idents))

CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat_fd@DB <- CellChatDB.use

cellchat_fd <- subsetData(cellchat_fd) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4)

cellchat_fd <- identifyOverExpressedGenes(cellchat_fd)
cellchat_fd <- identifyOverExpressedInteractions(cellchat_fd)  #LR

cellchat_fd <- computeCommunProb(cellchat_fd) #net
cellchat_fd <- filterCommunication(cellchat_fd, min.cells = 10)


cellchat_fd <- computeCommunProbPathway(cellchat_fd) #netP
cellchat_fd <- aggregateNet(cellchat_fd)

mat_fd <- cellchat_fd@net$weight

#subset fap only communication
#cellchat_fd@netP$prob[7:12,7:12,]<- 0
#cellchat_fd@netP$prob[7:12,7:12,]<- 0 
#cellchat_fd@net$pval[7:12,7:12,]<- 0
#cellchat_fd@net$count[7:12,7:12]<- 0 
#cellchat_fd@net$weight[7:12,7:12]<- 0

# Incoming/ Outgoing -----------------------------------------------------------
#PC
cellchat_pc<- netAnalysis_computeCentrality(cellchat_pc, slot.name = "netP")
selectK(cellchat_pc, pattern = "outgoing")
nPatterns <- 3 
cellchat_pc <- identifyCommunicationPatterns(cellchat_pc, pattern = "outgoing", k = nPatterns)

selectK(cellchat_pc, pattern = "incoming")
nPatterns <- 4
cellchat_pc <- identifyCommunicationPatterns(cellchat_pc, pattern = "incoming", k = nPatterns)

#PD
cellchat_pd<- netAnalysis_computeCentrality(cellchat_pd, slot.name = "netP")
selectK(cellchat_pd, pattern = "outgoing")
nPatterns <- 2 
cellchat_pd <- identifyCommunicationPatterns(cellchat_pd, pattern = "outgoing", k = nPatterns)

selectK(cellchat_pd, pattern = "incoming")
nPatterns <- 3
cellchat_pd <- identifyCommunicationPatterns(cellchat_pd, pattern = "incoming", k = nPatterns)

#FC
cellchat_fc<- netAnalysis_computeCentrality(cellchat_fc, slot.name = "netP")
selectK(cellchat_fc, pattern = "outgoing")
nPatterns <- 4 
cellchat_fc <- identifyCommunicationPatterns(cellchat_fc, pattern = "outgoing", k = nPatterns)

selectK(cellchat_fc, pattern = "incoming")
nPatterns <- 4
cellchat_fc <- identifyCommunicationPatterns(cellchat_fc, pattern = "incoming", k = nPatterns)

#FD
cellchat_fd<- netAnalysis_computeCentrality(cellchat_fd, slot.name = "netP")
selectK(cellchat_fd, pattern = "outgoing")
nPatterns <- 3 
cellchat_fd <- identifyCommunicationPatterns(cellchat_fd, pattern = "outgoing", k = nPatterns)

selectK(cellchat_fd, pattern = "incoming")
nPatterns <- 5
cellchat_fd <- identifyCommunicationPatterns(cellchat_fd, pattern = "incoming", k = nPatterns)

# Save Files --------------------------------------------------------------
saveRDS(cellchat_pc, "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Objects/specific_cellchat_pc.rds")
saveRDS(cellchat_pd, "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Objects/specific_cellchat_pd.rds")
saveRDS(cellchat_fc, "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Objects/specific_cellchat_fc.rds")
saveRDS(cellchat_fd, "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Objects/specific_cellchat_fd.rds")

# Load Files --------------------------------------------------------------
cellchat_pc <- readRDS("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Objects/specific_cellchat_pc.rds")
cellchat_pd <- readRDS("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Objects/specific_cellchat_pd.rds")
cellchat_fc <- readRDS("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Objects/specific_cellchat_fc.rds")
cellchat_fd <- readRDS("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Objects/specific_cellchat_fd.rds")

updateCellChat(cellchat_pc)
updateCellChat(cellchat_fc)
updateCellChat(cellchat_pd)
updateCellChat(cellchat_fd)

object.list <- list(`partial cuff` = cellchat_pc, 
                    `partial delt` = cellchat_pd,
                    `full cuff` = cellchat_fc, 
                    `full delt` = cellchat_fd)
cellchat <- mergeCellChat(object.list,add.names=names(object.list))
# Circle Plots ------------------------------------------------------------
pathways <- c("TGFb", "APP", "THY1", "VISFATIN", "TENASCIN", "SEMA3", "HGF", "BMP")
pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/HGF.pdf")

pathways.show <- pathways[7]
par(mfrow = c(2,2), xpd=TRUE)

weight.max <- getMaxWeight(object.list[c(1,3,4)], slot.name = c("netP"), attribute = pathways.show)
for (i in 1:length(object.list)) {
  tryCatch({
  netVisual_aggregate(object.list[[i]], vertex.receiver = vertex.receiver, signaling = pathways.show, layout = "circle", vertex.weight = NULL, edge.weight.max = weight.max[1], signaling.name = paste(pathways.show, names(object.list)[i]), arrow.size = 0.3)
  }, error=function(e){})
    }
dev.off()

#For when signaling is not showing up 
netVisual_aggregate(cellchat_pc, signaling = pathways.show, signaling.name = paste(pathways.show, edge.weight.max = weight.max[1],names(object.list)[1]),layout = "circle")

netVisual_aggregate(cellchat_pd, signaling = pathways.show,signaling.name = paste(pathways.show, edge.weight.max = weight.max[2], names(object.list)[2]), layout = "circle")

netVisual_aggregate(cellchat_fc, signaling = pathways.show, signaling.name = paste(pathways.show, edge.weight.max = weight.max[3], names(object.list)[3]),layout = "circle")

netVisual_aggregate(cellchat_fd, signaling = pathways.show,signaling.name = paste(pathways.show, edge.weight.max = weight.max[4], names(object.list)[4]), layout = "circle")

# LR Plots ---------------------------------------------------------
pathways <- c("TGFb", "APP", "THY1", "TENASCIN", "VISFATIN", "SEMA3", "HGF", "BMP")
pathways.show <- pathways

path <- "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/LR/"
for(i in 1:8){
  tryCatch({
  pdf(paste0(path,"specific_LR_",pathways.show[i],".pdf"))
  print(netAnalysis_contribution(cellchat_pc,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Partial Cuff"), signaling = pathways.show[i]))
  print(netAnalysis_contribution(cellchat_fc,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Full Cuff"), signaling = pathways.show[i]))
  print(netAnalysis_contribution(cellchat_pd,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Partial Delt"), signaling = pathways.show[i]))
  print(netAnalysis_contribution(cellchat_fd,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Full Delt"), signaling = pathways.show[i]))
  dev.off()
  }, error=function(e){})
}

#In case some objects do not have significant communication
pathways.show <- "HGF"
pdf(paste0("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/LR/specific_LR_",pathways.show,".pdf"))
netAnalysis_contribution(cellchat_pc,title = paste("Contribution of each L-R pair in", pathways.show, "Signaling Pathway, Partial Cuff"), signaling = pathways.show)
netAnalysis_contribution(cellchat_fc,title = paste("Contribution of each L-R pair in", pathways.show, "Signaling Pathway, Full Cuff"), signaling = pathways.show)
netAnalysis_contribution(cellchat_pd,title = paste("Contribution of each L-R pair in", pathways.show, "Signaling Pathway, Partial Delt"), signaling = pathways.show)
netAnalysis_contribution(cellchat_fd,title = paste("Contribution of each L-R pair in", pathways.show, "Signaling Pathway, Full Delt"), signaling = pathways.show)
dev.off()


# HeatMap -----------------------------------------------------------------
cellchat_pc <- netAnalysis_computeCentrality(cellchat_pc, slot.name = "netP") 
cellchat_pd <- netAnalysis_computeCentrality(cellchat_pd, slot.name = "netP")
cellchat_fc <- netAnalysis_computeCentrality(cellchat_fc, slot.name = "netP")
cellchat_fd <- netAnalysis_computeCentrality(cellchat_fd, slot.name = "netP")

pathways.show <- pathways

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Heatmap/specific_heatmap_pc.pdf")

#pc
for(i in 1:8){
  netAnalysis_signalingRole_network(cellchat_pc, 
                                    signaling = pathways.show[i], 
                                    width = 8, 
                                    height = 2.5, 
                                    font.size = 10)
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Heatmap/specific_heatmap_fc.pdf")
#fc
for(i in 1:8){
netAnalysis_signalingRole_network(cellchat_fc, 
                                  signaling = pathways.show[i], 
                                  width = 8, 
                                  height = 2.5, 
                                  font.size = 10)
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Heatmap/specific_heatmap_pd.pdf")
#pd
for(i in 1:8){
  tryCatch({
  netAnalysis_signalingRole_network(cellchat_pd, 
                                    signaling = pathways.show[i], 
                                    width = 8, 
                                    height = 2.5, 
                                    font.size = 10)
  }, error=function(e){})
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Heatmap/specific_heatmap_fd.pdf")
#fd
for(i in 1:8){
  netAnalysis_signalingRole_network(cellchat_fd, 
                                    signaling = pathways.show[i], 
                                    width = 8, 
                                    height = 2.5, 
                                    font.size = 10)
}
dev.off()


# NetVisual_bubble --------------------------------------------------------
netVisual_bubble <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, sort.by.source.priority = TRUE, color.heatmap = c("Spectral","viridis"), n.colors = 10, direction = -1, thresh = 0.05,
                             comparison = NULL, group = NULL, remove.isolate = FALSE, max.dataset = NULL, min.dataset = NULL,
                             min.quantile = 0, max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, color.text = NULL,
                             title.name = NULL, font.size = 10, font.size.title = 10, show.legend = TRUE,
                             grid.on = TRUE, color.grid = "grey90", angle.x = 90, vjust.x = NULL, hjust.x = NULL,
                             return.data = FALSE){
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  } else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 1, 1)
    vjust=c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
    })
  } else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      pairLR.use$pathway_name <- as.character(pairLR.use$pathway_name)
    } else if ("interaction_name" %in% colnames(pairLR.use)) {
      pairLR.use$interaction_name <- as.character(pairLR.use$interaction_name)
    }
  }  
  
  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }
    
    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)
    
    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T)*1.1, max(df.net$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    # rownames(df.net) <- df.net$interaction_name_2
    
    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")
    
    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2),])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, levels = cells.order)
    df <- df.net
  } else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net",
                                      sources.use = sources.use, targets.use = targets.use,
                                      signaling = signaling,
                                      pairLR.use = pairLR.use,
                                      thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      
      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], ")")
      
      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      } else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), ncol = 5))
        colnames(df.net) <- c("interaction_name_2","source.target","prob","pval","prob.original")
        df.net$source.target <- group.names0
      }
      # df.net$group.names <- sub(paste0(' \\(',dataset.name[comparison[ii]],'\\)'),'',as.character(df.net$source.target))
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T)*1.1, max(df.all$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    
    df <- df.all
    df <- with(df, df[order(interaction_name_2),])
    df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
    
    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i], " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }
  
  min.cutoff <- quantile(df$prob, min.quantile,na.rm= T)
  max.cutoff <- quantile(df$prob, max.quantile,na.rm= T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  
  
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    # line.on <- FALSE
    # df <- df[!is.na(df$prob),]
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, ,drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        #idx.na <- c(which(is.na(values)), which(!(dataset.name[comparison] %in% df.i.j$dataset)))
        dataset.na <- c(df.i.j$dataset[is.na(values)], setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          } else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            } else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
    #df <- df[!is.na(df$prob), ]
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  # Re-order y-axis
  if (!is.null(pairLR.use)) {
    interaction_name_2.order <- intersect(object@DB$interaction[pairLR.use$interaction_name, ]$interaction_name_2, unique(df$interaction_name_2))
    df$interaction_name_2 <- factor(df$interaction_name_2, levels = interaction_name_2.order)
  }
  
  # Re-order x-axis
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target),unique(df$source.target)))
  if (sort.by.target & !sort.by.source) {
    if (!is.null(targets.use)) {
      df$target <- factor(df$target, levels = intersect(targets.use, df$target))
      df <- with(df, df[order(target, source),])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & !sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, df$source))
      df <- with(df, df[order(source, target),])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, df$source))
      if (!is.null(targets.use)) {
        df$target <- factor(df$target, levels = intersect(targets.use, df$target))
      }
      if (sort.by.source.priority) {
        df <- with(df, df[order(source, target),])
      } else {
        df <- with(df, df[order(target, source),])
      }
      
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
    geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")
  
  values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
  g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
  #g <- g + scale_radius(range = c(1,3), breaks = values,labels = names(values), name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  } else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }
  
  g <- g + theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  
  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept=seq(1.5, length(unique(df$source.target))-0.5, 1),lwd=0.1,colour=color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour=color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5+length(dataset.name[comparison]), length(group.names0)*length(dataset.name[comparison]), by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype="dashed", color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      } else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      #names(color) <- dataset.name[comparison]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order)-1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  } else {
    return(g)
  }
  
}




# Fold Change -------------------------------------------------------------
pathways.show <- c("HGF","THY1","SEMA3", "BMP")
#PC vs FC
gg1<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(1,3),
                       title.name = "Increased signaling in Full Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 3,
                       return.data = T)


gg2<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(1,3),
                       title.name = "Decreased signaling in Full Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 1,
                       return.data = T)

write.csv(gg1$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/FC/specific_Inc_FCvPC.csv")
write.csv(gg2$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/FC/specific_Dec_FCvPC.csv")

#dot plot
pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/UpDown/fc_pc.pdf")

pathways <- pathways.show[4]


netVisual_bubble(cellchat, 
                 sources.use = c(1:12), 
                 targets.use = c(1:12),
                 pairLR.use = pairLR.use.up,
                 comparison = c(1,3),
                 title.name = paste0("Changed ",pathways, " signaling in Full Cuff vs Partial Cuff"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 3) + 
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0)

netVisual_bubble(cellchat, 
                 sources.use = c(1:12), 
                 targets.use = c(1:12),  
                 pairLR.use = pairLR.use.down,
                 comparison = c(1,3),
                 title.name = paste0("Decreased ",pathways, " signaling in Full Cuff vs Partial Cuff"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 1, 
                 return.data = T) +
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0.18) 

dev.off()


#PC vs PD
gg3<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(1,2),
                       title.name = "Increased signaling in Partial Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 1,
                       return.data = T)


gg4<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(1,2),
                       title.name = "Decreased signaling in Partial Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 2,
                       return.data = T)

write.csv(gg3$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/FC/specific_Inc_PCvPD.csv")
write.csv(gg4$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/FC/specific_Dec_PCvPD.csv")

#FC vs FD
gg5<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(3,4),
                       title.name = "Increased signaling in Full Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 3,
                       return.data = T)

gg6<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(3,4),
                       title.name = "Decreased signaling in Full Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 4,
                       return.data = T)
write.csv(gg5$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/FC/specific_Inc_FCvFD.csv")
write.csv(gg6$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/FC/specific_Dec_FCvFD.csv")


#dot plot
pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/UpDown/fc_fd.pdf")
  
  pathways <- pathways.show[4]
  
  netVisual_bubble(cellchat, 
                   sources.use = c(1:12), 
                   targets.use = c(1:12),
                   signaling = pathways,
                   comparison = c(3,4),
                   title.name = paste0("Increased ",pathways, " signaling in Full Cuff vs Full Delt"), 
                   angle.x = 45, 
                   remove.isolate = T,
                   max.dataset = 3) +
    scale_colour_gradient2(high = "#8b0000",                                  
                           low = "#B6D0E2",
                           midpoint = 0.162)
  
  netVisual_bubble(cellchat, 
                   sources.use = c(1:12), 
                   targets.use = c(1:12),  
                   signaling = pathways,
                   comparison = c(3,4),
                   title.name = paste0("Decreased ",pathways, " signaling in Full Cuff vs Full Delt"), 
                   angle.x = 45, 
                   remove.isolate = T,
                   max.dataset = 4) +
    scale_colour_gradient2(high = "#8b0000",                                  
                           low = "#B6D0E2",
                           midpoint = 0.19)
dev.off()

# Chord -------------------------------------------------------------------


pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Chord/specific_chord.pdf")
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]]@net$count, net = object.list[[i]]@net$count, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/Chord/specific_chord_weight.pdf")
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]]@net$weight, net = object.list[[i]]@net$weight, title.name = paste0("Interaction weights/strength - ", names(object.list)[i]))
}
dev.off()
