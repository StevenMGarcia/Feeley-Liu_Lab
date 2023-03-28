
# libraries ---------------------------------------------------------------
library(Seurat)
library(CellChat)
library(patchwork)
library(plyr)
library(NMF)
library(ggalluvial)
library(tidyr)
library(stringr)



# Load/ Combine Objects ---------------------------------------------------

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




# FC v PC: BMP ------------------------------------------------------------
pos.dataset = "full cuff"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "full cuff",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL, 
                              signaling = "BMP")
net.down <- subsetCommunication(cellchat, 
                                net = net, 
                                datasets = "partial cuff",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1,
                                signaling = "BMP")

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pathways.show <- c("HGF","THY1","SEMA3", "BMP")
pathways <- pathways.show[4]

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c(1:12), 
                        targets.use = c(1:12),
                        pairLR.use = pairLR.use.up,
                        comparison = c(1,3),
                        title.name = paste0("Increased ",pathways, " signaling in Full Cuff Compared to Partial Cuff"), 
                        angle.x = 45, 
                        remove.isolate = T,
                        max.dataset = 3) + 
  scale_colour_gradient2(high = "#8b0000",
                         midpoint = 0)
BMP <- ggplot_build(gg1)

df_bmp <- BMP$plot$data


df_bmp_fc <- subset(df_bmp, dataset == "full cuff")
df_bmp_fc$prob_fc <- df_bmp_fc$prob
df_bmp_fd <- subset(df_bmp, dataset == "partial cuff")
df_bmp_fd$prob_pc <- df_bmp_fd$prob

df_bmp <- rbind.fill(df_bmp_fd, df_bmp_fc)
df_bmp$source1 <- df_bmp$source
df_bmp$target1 <- df_bmp$target
df_bmp$ligand1 <- df_bmp$ligand
df_bmp$receptor1 <- df_bmp$receptor
df_bmp <- unite(df_bmp, col='combo', c('source1', 'target1',"ligand1", "receptor1"), sep='')

df_bmp<- df_bmp %>%
  group_by(combo) %>% summarise_each(funs(first(.[!is.na(.)])))

df_bmp[is.na(df_bmp)] <- 0



df_bmp$prob <- df_bmp$prob_fc - df_bmp$prob_pc
df_bmp$source.target <- gsub("\\s*\\([^\\)]+\\)","",as.character(df_bmp$source.target))

df_bmp$combo <- NULL
df_bmp$prob_fc <- NULL
df_bmp$prob_pc <- NULL

BMP$plot$data <- df_bmp

BMP$plot$data$dataset <- NULL

BMP_fcpc <- BMP
BMP_fcpc

# FC v FD: BMP  -----------------------------------------------------------
pos.dataset = "full cuff"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "full cuff",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL, 
                              signaling = "BMP")
net.down <- subsetCommunication(cellchat, 
                                net = net, 
                                datasets = "full delt",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1,
                                signaling = "BMP")

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pathways.show <- c("HGF","THY1","SEMA3", "BMP")
pathways <- pathways.show[4]

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c(1:12), 
                        targets.use = c(1:12),
                        pairLR.use = pairLR.use.up,
                        comparison = c(3,4),
                        title.name = paste0("Increased ",pathways, " signaling in Full Cuff Compared to Full Delt"), 
                        angle.x = 45, 
                        remove.isolate = T,
                        max.dataset = 3,
                        return.data = T) + 
  scale_colour_gradient2(high = "#8b0000",
                         midpoint = 0.00)
BMP <- ggplot_build(gg1)

df_bmp <- BMP$plot$data


df_bmp_fc <- subset(df_bmp, dataset == "full cuff")
df_bmp_fc$prob_fc <- df_bmp_fc$prob
df_bmp_fd <- subset(df_bmp, dataset == "full delt")
df_bmp_fd$prob_fd <- df_bmp_fd$prob

df_bmp <- rbind.fill(df_bmp_fd, df_bmp_fc)
df_bmp$source1 <- df_bmp$source
df_bmp$target1 <- df_bmp$target
df_bmp$ligand1 <- df_bmp$ligand
df_bmp$receptor1 <- df_bmp$receptor
df_bmp <- unite(df_bmp, col='combo', c('source1', 'target1',"ligand1", "receptor1"), sep='')

df_bmp<- df_bmp %>%
  group_by(combo) %>% summarise_each(funs(first(.[!is.na(.)])))

df_bmp[is.na(df_bmp)] <- 0



df_bmp$prob <- df_bmp$prob_fc - df_bmp$prob_fd
df_bmp$source.target <- gsub("\\s*\\([^\\)]+\\)","",as.character(df_bmp$source.target))

df_bmp$combo <- NULL
df_bmp$prob_fc <- NULL
df_bmp$prob_fd <- NULL

BMP$plot$data <- df_bmp

BMP$plot$data$dataset <- NULL

BMP_fcfd <- BMP

BMP_fcfd





# FC v PC: SEMA3 ------------------------------------------------------------
pos.dataset = "full cuff"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "full cuff",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL, 
                              signaling = "SEMA3")
net.down <- subsetCommunication(cellchat, 
                                net = net, 
                                datasets = "partial cuff",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1,
                                signaling = "SEMA3")

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pathways.show <- c("HGF","THY1","SEMA3", "BMP")
pathways <- pathways.show[3]

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c(1:12), 
                        targets.use = c(1:12),
                        pairLR.use = pairLR.use.up,
                        comparison = c(1,3),
                        title.name = paste0("Increased ",pathways, " signaling in Full Cuff Compared to Partial Cuff"), 
                        angle.x = 45, 
                        remove.isolate = T,
                        max.dataset = 3) + 
  scale_colour_gradient2(high = "#8b0000",
                         midpoint = 0)
SEMA3 <- ggplot_build(gg1)

df_sema3 <- SEMA3$plot$data


df_sema3_fc <- subset(df_sema3, dataset == "full cuff")
df_sema3_fc$prob_fc <- df_sema3_fc$prob
df_sema3_fd <- subset(df_sema3, dataset == "partial cuff")
df_sema3_fd$prob_pc <- df_sema3_fd$prob

df_sema3 <- rbind.fill(df_sema3_fd, df_sema3_fc)
df_sema3$source1 <- df_sema3$source
df_sema3$target1 <- df_sema3$target
df_sema3$ligand1 <- df_sema3$ligand
df_sema3$receptor1 <- df_sema3$receptor
df_sema3 <- unite(df_sema3, col='combo', c('source1', 'target1',"ligand1", "receptor1"), sep='')

df_sema3<- df_sema3 %>%
  group_by(combo) %>% summarise_each(funs(first(.[!is.na(.)])))

df_sema3[is.na(df_sema3)] <- 0



df_sema3$prob <- df_sema3$prob_fc - df_sema3$prob_pc
df_sema3$source.target <- gsub("\\s*\\([^\\)]+\\)","",as.character(df_sema3$source.target))

df_sema3$combo <- NULL
df_sema3$prob_fc <- NULL
df_sema3$prob_pc <- NULL

SEMA3$plot$data <- df_sema3

SEMA3$plot$data$dataset <- NULL

SEMA3_fcpc <- SEMA3
SEMA3_fcpc

# FC v FD: SEMA3  -----------------------------------------------------------
pos.dataset = "full cuff"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "full cuff",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL, 
                              signaling = "SEMA3")
net.down <- subsetCommunication(cellchat, 
                                net = net, 
                                datasets = "full delt",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1,
                                signaling = "SEMA3")

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pathways.show <- c("HGF","THY1","SEMA3", "BMP")
pathways <- pathways.show[3]

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c(1:12), 
                        targets.use = c(1:12),
                        pairLR.use = pairLR.use.up,
                        comparison = c(3,4),
                        title.name = paste0("Increased ",pathways, " signaling in Full Cuff Compared to Full Delt"), 
                        angle.x = 45, 
                        remove.isolate = T,
                        max.dataset = 3) + 
  scale_colour_gradient2(high = "#8b0000",
                         midpoint = 0.00)
SEMA3 <- ggplot_build(gg1)

df_sema3 <- SEMA3$plot$data


df_sema3_fc <- subset(df_sema3, dataset == "full cuff")
df_sema3_fc$prob_fc <- df_sema3_fc$prob
df_sema3_fd <- subset(df_sema3, dataset == "full delt")
df_sema3_fd$prob_fd <- df_sema3_fd$prob

df_sema3 <- rbind.fill(df_sema3_fd, df_sema3_fc)
df_sema3$source1 <- df_sema3$source
df_sema3$target1 <- df_sema3$target
df_sema3$ligand1 <- df_sema3$ligand
df_sema3$receptor1 <- df_sema3$receptor
df_sema3 <- unite(df_sema3, col='combo', c('source1', 'target1',"ligand1", "receptor1"), sep='')

df_sema3<- df_sema3 %>%
  group_by(combo) %>% summarise_each(funs(first(.[!is.na(.)])))

df_sema3[is.na(df_sema3)] <- 0



df_sema3$prob <- df_sema3$prob_fc - df_sema3$prob_fd
df_sema3$source.target <- gsub("\\s*\\([^\\)]+\\)","",as.character(df_sema3$source.target))

df_sema3$combo <- NULL
df_sema3$prob_fc <- NULL
df_sema3$prob_fd <- NULL

SEMA3$plot$data <- df_sema3

SEMA3$plot$data$dataset <- NULL

SEMA3_fcfd <- SEMA3
                        
SEMA3_fcfd
                        
                        
                        





# FC v PC: THY1 ------------------------------------------------------------
pos.dataset = "full cuff"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "full cuff",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL, 
                              signaling = "THY1")
net.down <- subsetCommunication(cellchat, 
                                net = net, 
                                datasets = "partial cuff",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1,
                                signaling = "THY1")

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pathways.show <- c("HGF","THY1","SEMA3", "BMP")
pathways <- pathways.show[2]

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c(1:12), 
                        targets.use = c(1:12),
                        pairLR.use = pairLR.use.up,
                        comparison = c(1,3),
                        title.name = paste0("Increased ",pathways, " signaling in Full Cuff Compared to Partial Cuff"), 
                        angle.x = 45, 
                        remove.isolate = T,
                        max.dataset = 3) + 
  scale_colour_gradient2(high = "#8b0000",
                         midpoint = 0)
THY1 <- ggplot_build(gg1)

df_thy1 <- THY1$plot$data


df_thy1_fc <- subset(df_thy1, dataset == "full cuff")
df_thy1_fc$prob_fc <- df_thy1_fc$prob
df_thy1_fd <- subset(df_thy1, dataset == "partial cuff")
df_thy1_fd$prob_pc <- df_thy1_fd$prob

df_thy1 <- rbind.fill(df_thy1_fd, df_thy1_fc)
df_thy1$source1 <- df_thy1$source
df_thy1$target1 <- df_thy1$target
df_thy1$ligand1 <- df_thy1$ligand
df_thy1$receptor1 <- df_thy1$receptor
df_thy1 <- unite(df_thy1, col='combo', c('source1', 'target1',"ligand1", "receptor1"), sep='')

df_thy1<- df_thy1 %>%
  group_by(combo) %>% summarise_each(funs(first(.[!is.na(.)])))

df_thy1[is.na(df_thy1)] <- 0



df_thy1$prob <- df_thy1$prob_fc - df_thy1$prob_pc
df_thy1$source.target <- gsub("\\s*\\([^\\)]+\\)","",as.character(df_thy1$source.target))

df_thy1$combo <- NULL
df_thy1$prob_fc <- NULL
df_thy1$prob_pc <- NULL

THY1$plot$data <- df_thy1

THY1$plot$data$dataset <- NULL

THY1_fcpc <- THY1
THY1_fcpc

# FC v FD: THY1  -----------------------------------------------------------
pos.dataset = "full cuff"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "full cuff",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL, 
                              signaling = "THY1")
net.down <- subsetCommunication(cellchat, 
                                net = net, 
                                datasets = "full delt",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1,
                                signaling = "THY1")

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pathways.show <- c("HGF","THY1","SEMA3", "BMP")
pathways <- pathways.show[2]

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c(1:12), 
                        targets.use = c(1:12),
                        pairLR.use = pairLR.use.up,
                        comparison = c(3,4),
                        title.name = paste0("Increased ",pathways, " signaling in Full Cuff Compared to Full Delt"), 
                        angle.x = 45, 
                        remove.isolate = T,
                        max.dataset = 3) + 
  scale_colour_gradient2(high = "#8b0000",
                         midpoint = 0.00)
THY1 <- ggplot_build(gg1)

df_thy1 <- THY1$plot$data


df_thy1_fc <- subset(df_thy1, dataset == "full cuff")
df_thy1_fc$prob_fc <- df_thy1_fc$prob
df_thy1_fd <- subset(df_thy1, dataset == "full delt")
df_thy1_fd$prob_fd <- df_thy1_fd$prob

df_thy1 <- rbind.fill(df_thy1_fd, df_thy1_fc)
df_thy1$source1 <- df_thy1$source
df_thy1$target1 <- df_thy1$target
df_thy1$ligand1 <- df_thy1$ligand
df_thy1$receptor1 <- df_thy1$receptor
df_thy1 <- unite(df_thy1, col='combo', c('source1', 'target1',"ligand1", "receptor1"), sep='')

df_thy1<- df_thy1 %>%
  group_by(combo) %>% summarise_each(funs(first(.[!is.na(.)])))

df_thy1[is.na(df_thy1)] <- 0



df_thy1$prob <- df_thy1$prob_fc - df_thy1$prob_fd
df_thy1$source.target <- gsub("\\s*\\([^\\)]+\\)","",as.character(df_thy1$source.target))

df_thy1$combo <- NULL
df_thy1$prob_fc <- NULL
df_thy1$prob_fd <- NULL

THY1$plot$data <- df_thy1

THY1$plot$data$dataset <- NULL

THY1_fcfd <- THY1

THY1_fcfd


# FC v PC: HGF ------------------------------------------------------------
pos.dataset = "full cuff"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "full cuff",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL, 
                              signaling = "HGF")
net.down <- subsetCommunication(cellchat, 
                                net = net, 
                                datasets = "partial cuff",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1,
                                signaling = "HGF")

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pathways.show <- c("HGF","THY1","SEMA3", "BMP")
pathways <- pathways.show[1]

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c(1:12), 
                        targets.use = c(1:12),
                        pairLR.use = pairLR.use.up,
                        comparison = c(1,3),
                        title.name = paste0("Increased ",pathways, " signaling in Full Cuff Compared to Partial Cuff"), 
                        angle.x = 45, 
                        remove.isolate = T,
                        max.dataset = 3) + 
  scale_colour_gradient2(high = "#8b0000",
                         midpoint = 0)
HGF <- ggplot_build(gg1)

df_hgf <- HGF$plot$data


df_hgf_fc <- subset(df_hgf, dataset == "full cuff")
df_hgf_fc$prob_fc <- df_hgf_fc$prob
df_hgf_fd <- subset(df_hgf, dataset == "partial cuff")
df_hgf_fd$prob_pc <- df_hgf_fd$prob

df_hgf <- rbind.fill(df_hgf_fd, df_hgf_fc)
df_hgf$source1 <- df_hgf$source
df_hgf$target1 <- df_hgf$target
df_hgf$ligand1 <- df_hgf$ligand
df_hgf$receptor1 <- df_hgf$receptor
df_hgf <- unite(df_hgf, col='combo', c('source1', 'target1',"ligand1", "receptor1"), sep='')

df_hgf<- df_hgf %>%
  group_by(combo) %>% summarise_each(funs(first(.[!is.na(.)])))

df_hgf[is.na(df_hgf)] <- 0



df_hgf$prob <- df_hgf$prob_fc - df_hgf$prob_pc
df_hgf$source.target <- gsub("\\s*\\([^\\)]+\\)","",as.character(df_hgf$source.target))

df_hgf$combo <- NULL
df_hgf$prob_fc <- NULL
df_hgf$prob_pc <- NULL

HGF$plot$data <- df_hgf

HGF$plot$data$dataset <- NULL

HGF_fcpc <- HGF
HGF_fcpc

# FC v FD: HGF  -----------------------------------------------------------
pos.dataset = "full cuff"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "full cuff",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL, 
                              signaling = "HGF")
net.down <- subsetCommunication(cellchat, 
                                net = net, 
                                datasets = "full delt",
                                ligand.logFC = -0.1, 
                                receptor.logFC = -0.1,
                                signaling = "HGF")

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pathways.show <- c("HGF","THY1","SEMA3", "BMP")
pathways <- pathways.show[1]

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c(1:12), 
                        targets.use = c(1:12),
                        pairLR.use = pairLR.use.up,
                        comparison = c(3,4),
                        title.name = paste0("Increased ",pathways, " signaling in Full Cuff Compared to Full Delt"), 
                        angle.x = 45, 
                        remove.isolate = T,
                        max.dataset = 3) + 
  scale_colour_gradient2(high = "#8b0000",
                         midpoint = 0.00)
HGF <- ggplot_build(gg1)

df_hgf <- HGF$plot$data


df_hgf_fc <- subset(df_hgf, dataset == "full cuff")
df_hgf_fc$prob_fc <- df_hgf_fc$prob
df_hgf_fd <- subset(df_hgf, dataset == "full delt")
df_hgf_fd$prob_fd <- df_hgf_fd$prob

df_hgf <- rbind.fill(df_hgf_fd, df_hgf_fc)
df_hgf$source1 <- df_hgf$source
df_hgf$target1 <- df_hgf$target
df_hgf$ligand1 <- df_hgf$ligand
df_hgf$receptor1 <- df_hgf$receptor
df_hgf <- unite(df_hgf, col='combo', c('source1', 'target1',"ligand1", "receptor1"), sep='')

df_hgf<- df_hgf %>%
  group_by(combo) %>% summarise_each(funs(first(.[!is.na(.)])))

df_hgf[is.na(df_hgf)] <- 0



df_hgf$prob <- df_hgf$prob_fc - df_hgf$prob_fd
df_hgf$source.target <- gsub("\\s*\\([^\\)]+\\)","",as.character(df_hgf$source.target))

df_hgf$combo <- NULL
df_hgf$prob_fc <- NULL
df_hgf$prob_fd <- NULL

HGF$plot$data <- df_hgf

HGF$plot$data$dataset <- NULL

HGF_fcfd <- HGF

HGF_fcfd


# Axis text ---------------------------------------------------------------
BMP_fcpc[["plot"]][["theme"]][["axis.text.x"]][["colour"]] <- "black"
BMP_fcfd[["plot"]][["theme"]][["axis.text.x"]][["colour"]] <- "black"
  
SEMA3_fcpc[["plot"]][["theme"]][["axis.text.x"]][["colour"]] <- "black"
SEMA3_fcfd[["plot"]][["theme"]][["axis.text.x"]][["colour"]] <- "black"

THY1_fcpc[["plot"]][["theme"]][["axis.text.x"]][["colour"]] <- "black"
THY1_fcfd[["plot"]][["theme"]][["axis.text.x"]][["colour"]] <- "black"

HGF_fcpc[["plot"]][["theme"]][["axis.text.x"]][["colour"]] <- "black"
HGF_fcfd[["plot"]][["theme"]][["axis.text.x"]][["colour"]] <- "black"


# Heatmap conversion ------------------------------------------------------
ggplot(data = BMP_fcfd$plot$data, aes(x = source.target, y = interaction_name_2)) +
  geom_tile(aes(fill = prob)) +
  scale_fill_gradient(high = "purple",
                       low = "yellow") +
  theme(axis.text.y = element_text(), 
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.7),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())+
  labs(fill = "Commun Prob")+
  ggtitle("Increased BMP signaling in Full Cuff Compared to Full Delt")+
  coord_equal()
  

# Plot --------------------------------------------------------------------



rep_str = c("GLI-1"='GLI1','Satellite'='Myogenic','SCX/ TNMD'='SCX/TNMD')
BMP_fcpc$plot$data$source.target <- str_replace_all(BMP_fcpc$plot$data$source.target, rep_str)
BMP_fcfd$plot$data$source.target <- str_replace_all(BMP_fcfd$plot$data$source.target, rep_str)
SEMA3_fcpc$plot$data$source.target <- str_replace_all(SEMA3_fcpc$plot$data$source.target, rep_str)
SEMA3_fcfd$plot$data$source.target <- str_replace_all(SEMA3_fcfd$plot$data$source.target, rep_str)


pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/Specific/UpDown/Heatmap.pdf",width=15, height=5)
ggplot(data = BMP_fcpc$plot$data, aes(x = source.target, y = interaction_name_2)) +
  geom_tile(aes(fill = prob)) +
  scale_fill_gradient(high = "yellow",
                      low = "purple") +
  theme(axis.text.y = element_text(), 
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.7),
        plot.title = element_text(hjust = 0.5, vjust = 3),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())+
  labs(fill = "Commun Prob")+
  ggtitle("Increased BMP signaling in Full Cuff Compared to Partial Cuff")+
  coord_equal()


ggplot(data = BMP_fcfd$plot$data, aes(x = source.target, y = interaction_name_2)) +
  geom_tile(aes(fill = prob)) +
  scale_fill_gradient(high = "yellow",
                      low = "purple") +
  theme(axis.text.y = element_text(), 
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.7),
        plot.title = element_text(hjust = 0.5, vjust = 3),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())+
  labs(fill = "Commun Prob")+
  ggtitle("Increased BMP signaling in Full Cuff Compared to Full Delt")+
  coord_equal()

ggplot(data = SEMA3_fcpc$plot$data, aes(x = source.target, y = interaction_name_2)) +
  geom_tile(aes(fill = prob)) +
  scale_fill_gradient(high = "yellow",
                      low = "purple") +
  theme(axis.text.y = element_text(), 
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.7),
        plot.title = element_text(hjust = 0.5, vjust = 3),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())+
  labs(fill = "Commun Prob")+
  ggtitle("Increased SEMA3 signaling in Full Cuff Compared to Partial Cuff")+
  coord_equal()


ggplot(data = SEMA3_fcfd$plot$data, aes(x = source.target, y = interaction_name_2)) +
  geom_tile(aes(fill = prob)) +
  scale_fill_gradient(high = "yellow",
                      low = "purple") +
  theme(axis.text.y = element_text(), 
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.7),
        plot.title = element_text(hjust = 0.5, vjust = 3),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())+
  labs(fill = "Commun Prob")+
  ggtitle("Increased SEMA3 signaling in Full Cuff Compared to Full Delt")+
  coord_equal()

dev.off()

THY1_fcpc
THY1_fcfd

HGF_fcpc
HGF_fcfd

dev.off()

