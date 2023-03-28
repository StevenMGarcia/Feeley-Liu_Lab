# Libraries ---------------------------------------------------------------
library(Seurat)
library(CellChat)
library(patchwork)
library(plyr)
library(NMF)
library(ggalluvial)
library(devtools)


# Load Object ------------------------------------------------------
integrated <- readRDS("~/Desktop/UCSF/Steven/CellChat/hum_muscle_hetero_faps.rds") 

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


levels(integrated@meta.data$reanalyzed) <-c("FAP","FAP","FAP","FAP","FAP","FAP","FAP","FAP","FAP","FAP","FAP", "Fibrogenic","Pericytes","Endothelial", "Satellite", "Macrophages", "T-Cells")
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

cellchat_pc <- computeCommunProb(cellchat_pc) #net
cellchat_pc <- filterCommunication(cellchat_pc, min.cells = 10)


cellchat_pc <- computeCommunProbPathway(cellchat_pc) #netP
cellchat_pc <- aggregateNet(cellchat_pc)


#subset fap only communication
#cellchat_pc@netP$prob
#cellchat_pc@netP$prob[2:7,2:7,]<- 0
#cellchat_pc@net$prob[2:7,2:7,] <- 0 
#cellchat_pc@net$pval[2:7,2:7,] <- 0
#cellchat_pc@net$count[2:7,2:7] <- 0 
#cellchat_pc@net$weight[2:7,2:7] <- 0


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

#subset fap only communication
cellchat_pd <- computeCommunProbPathway(cellchat_pd) #netP
cellchat_pd <- aggregateNet(cellchat_pd)

#cellchat_pd@netP$prob
#cellchat_pd@netP$prob[2:7,2:7,]<- 0
#cellchat_pd@net$prob[2:7,2:7,] <- 0 
#cellchat_pd@net$pval[2:7,2:7,] <- 0
#cellchat_pd@net$count[2:7,2:7] <- 0 
#cellchat_pd@net$weight[2:7,2:7] <- 0

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
cellchat_fc@netP$prob
cellchat_fc <- aggregateNet(cellchat_fc)


#subset fap only communication
#cellchat_fc@netP$prob[2:7,2:7,]<- 0
#cellchat_fc@net$prob[2:7,2:7,] <- 0 
#cellchat_fc@net$pval[2:7,2:7,] <- 0
#cellchat_fc@net$count[2:7,2:7] <- 0 
#cellchat_fc@net$weight[2:7,2:7] <- 0

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
cellchat_fd@netP$prob

cellchat_fd <- aggregateNet(cellchat_fd)


#subset fap only communication
#cellchat_fd@netP$prob[2:7,2:7,]<- 0
#cellchat_fd@net$prob[2:7,2:7,] <- 0 
#cellchat_fd@net$pval[2:7,2:7,] <- 0
#cellchat_fd@net$count[2:7,2:7] <- 0 
#cellchat_fd@net$weight[2:7,2:7] <- 0

# Incoming/ Outgoing -----------------------------------------------------------
#PC
cellchat_pc<- netAnalysis_computeCentrality(cellchat_pc, slot.name = "netP")
selectK(cellchat_pc, pattern = "outgoing")
nPatterns <- 3 
cellchat_pc <- identifyCommunicationPatterns(cellchat_pc, pattern = "outgoing", k = nPatterns)

selectK(cellchat_pc, pattern = "incoming")
nPatterns <- 5
cellchat_pc <- identifyCommunicationPatterns(cellchat_pc, pattern = "incoming", k = nPatterns)

#PD
cellchat_pd<- netAnalysis_computeCentrality(cellchat_pd, slot.name = "netP")
selectK(cellchat_pd, pattern = "outgoing")
nPatterns <- 4 
cellchat_pd <- identifyCommunicationPatterns(cellchat_pd, pattern = "outgoing", k = nPatterns)

selectK(cellchat_pd, pattern = "incoming")
nPatterns <- 4
cellchat_pd <- identifyCommunicationPatterns(cellchat_pd, pattern = "incoming", k = nPatterns)

#FC
cellchat_fc<- netAnalysis_computeCentrality(cellchat_fc, slot.name = "netP")
selectK(cellchat_fc, pattern = "outgoing")
nPatterns <- 3 
cellchat_fc <- identifyCommunicationPatterns(cellchat_fc, pattern = "outgoing", k = nPatterns)

selectK(cellchat_fc, pattern = "incoming")
nPatterns <- 4
cellchat_fc <- identifyCommunicationPatterns(cellchat_fc, pattern = "incoming", k = nPatterns)

#FD
cellchat_fd<- netAnalysis_computeCentrality(cellchat_fd, slot.name = "netP")
selectK(cellchat_fd, pattern = "outgoing")
nPatterns <- 4 
cellchat_fd <- identifyCommunicationPatterns(cellchat_fd, pattern = "outgoing", k = nPatterns)

selectK(cellchat_fd, pattern = "incoming")
nPatterns <- 3
cellchat_fd <- identifyCommunicationPatterns(cellchat_fd, pattern = "incoming", k = nPatterns)

# Save Files --------------------------------------------------------------
saveRDS(cellchat_pc, "~/Desktop/UCSF/Steven/Cellchat/Paper/General/Objects/gen_cellchat_pc.rds")
saveRDS(cellchat_pd, "~/Desktop/UCSF/Steven/Cellchat/Paper/General/Objects/gen_cellchat_pd.rds")
saveRDS(cellchat_fc, "~/Desktop/UCSF/Steven/Cellchat/Paper/General/Objects/gen_cellchat_fc.rds")
saveRDS(cellchat_fd, "~/Desktop/UCSF/Steven/Cellchat/Paper/General/Objects/gen_cellchat_fd.rds")

# Load Files ------------------------------------------------------------
cellchat_pc <- readRDS("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Objects/gen_cellchat_pc.rds")
cellchat_pd <- readRDS("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Objects/gen_cellchat_pd.rds")
cellchat_fc <- readRDS("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Objects/gen_cellchat_fc.rds")
cellchat_fd <- readRDS("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Objects/gen_cellchat_fd.rds")

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
pathways <- c("TGFb", "APP", "THY1", "TENASCIN", "VISFATIN")

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/general_circle.pdf")

pathways.show <- pathways[5]
par(mfrow = c(2,2), xpd=TRUE)

weight.max <- getMaxWeight(object.list[1:4], slot.name = c("netP"), attribute = pathways.show)
for (i in 1:length(object.list)) {
  tryCatch({
    netVisual_aggregate(object.list[[i]], vertex.receiver = vertex.receiver, signaling = pathways.show, layout = "circle", vertex.weight = NULL, edge.weight.max = weight.max[1], signaling.name = paste(pathways.show, names(object.list)[i]), arrow.size = 0.3)
  }, error=function(e){})
}

dev.off()

#When signaling for certain groups does not show up 
netVisual_aggregate(cellchat_pc, signaling = pathways.show, signaling.name = paste(pathways.show, names(object.list)[1]),layout = "circle")

netVisual_aggregate(cellchat_pd, signaling = pathways.show,signaling.name = paste(pathways.show, names(object.list)[2]), layout = "circle")

netVisual_aggregate(cellchat_fc, signaling = pathways.show, signaling.name = paste(pathways.show, names(object.list)[3]),layout = "circle")

netVisual_aggregate(cellchat_fd, signaling = pathways.show,signaling.name = paste(pathways.show, names(object.list)[4]), layout = "circle")


# LR Plots ----------------------------------------------------------------------
pathways <- c("TGFb", "APP", "THY1", "TENASCIN", "VISFATIN", "SEMA3", "HGF", "BMP")
pathways.show <- pathways

path <- "~/Desktop/UCSF/Steven/Cellchat/Paper/General/LR/"
for(i in 1:8){
  tryCatch({
    pdf(paste0(path,"general_LR_",pathways.show[i],".pdf"))
    print(netAnalysis_contribution(cellchat_pc,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Partial Cuff"), signaling = pathways.show[i]))
    print(netAnalysis_contribution(cellchat_fc,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Full Cuff"), signaling = pathways.show[i]))
    print(netAnalysis_contribution(cellchat_pd,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Partial Delt"), signaling = pathways.show[i]))
    print(netAnalysis_contribution(cellchat_fd,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Full Delt"), signaling = pathways.show[i]))
    dev.off()
  }, error=function(e){})
}

#In case some objects do not have significant communication
pathways.show <- "HGF"
pdf(paste0("~/Desktop/UCSF/Steven/Cellchat/Paper/General/LR/general_LR_",pathways.show,".pdf"))
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

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Heatmap/general_heatmap_pc.pdf")

#pc
for(i in 1:8){
  netAnalysis_signalingRole_network(cellchat_pc, 
                                    signaling = pathways.show[i], 
                                    width = 8, 
                                    height = 2.5, 
                                    font.size = 10)
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Heatmap/general_heatmap_fc.pdf")
#fc
for(i in 1:8){
  netAnalysis_signalingRole_network(cellchat_fc, 
                                    signaling = pathways.show[i], 
                                    width = 8, 
                                    height = 2.5, 
                                    font.size = 10)
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Heatmap/general_heatmap_pd.pdf")
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

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Heatmap/general_heatmap_fd.pdf")
#fd
for(i in 1:8){
  netAnalysis_signalingRole_network(cellchat_fd, 
                                    signaling = pathways.show[i], 
                                    width = 8, 
                                    height = 2.5, 
                                    font.size = 10)
}
dev.off()

# Fold Change -------------------------------------------------------------
pathways.show <- c("HGF", "THY1", "SEMA3","BMP")
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

write.csv(gg1$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Inc_FCvPC.csv")
write.csv(gg2$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Dec_FCvPC.csv")

#dot plot
pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/UpDown/fc_pc.pdf")
pathways <- pathways.show[4]



netVisual_bubble(cellchat, 
                 sources.use = c(1:7), 
                 targets.use = c(1:7),
                 signaling = pathways,
                 comparison = c(1,3),
                 color.heatmap = c("viridis"),
                 title.name = paste0("Increased ",pathways, " signaling in Full Cuff vs Partial Cuff"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 3) +
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0.16)

netVisual_bubble(cellchat, 
                 sources.use = c(1:7), 
                 targets.use = c(1:7),  
                 signaling = pathways,
                 comparison = c(1,3),
                 title.name = paste0("Decreased ",pathways, " signaling in Full Cuff vs Partial Cuff"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 1) +
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0.137)
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

write.csv(gg3$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Inc_PCvPD.csv")
write.csv(gg4$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Dec_PCvPD.csv")

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

write.csv(gg5$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Inc_FCvFD.csv")
write.csv(gg6$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Dec_FCvFD.csv")

#dot plot
pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/UpDown/fc_fd.pdf")

pathways <- pathways.show[1]

netVisual_bubble(cellchat, 
                 sources.use = c(1:7), 
                 targets.use = c(1:7),
                 signaling = pathways,
                 comparison = c(3,4),
                 title.name = paste0("Increased ",pathways, " signaling in Full Cuff vs Full Delt"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 3) +
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0.17)

netVisual_bubble(cellchat, 
                 sources.use = c(1:7), 
                 targets.use = c(1:7),  
                 signaling = pathways,
                 comparison = c(3,4),
                 title.name = paste0("Decreased ",pathways, " signaling in Full Cuff vs Full Delt"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 4) +
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0.17)
dev.off()



# Chord -------------------------------------------------------------------

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Chord/general_chord.pdf")
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]]@net$count, net = object.list[[i]]@net$count, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Chord/general_chord_weight.pdf")
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]]@net$weight, net = object.list[[i]]@net$weight, title.name = paste0("Interaction weights/strength - ", names(object.list)[i]))
}
dev.off()
