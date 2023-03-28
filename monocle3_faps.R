suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(dplyr)
  library(viridis)
  library(reshape2)
  library(batchelor)
  library(ggplot2)
  library(ggrastr)
  library(multiClust)
})

#### Create a Monocle CDS Object
# Project PC dimensions to whole data set
my.so=readRDS(file="../only_faps_integrated.rds")
Idents(my.so)<-"integrated_snn_res.0.3"
my.so <- ProjectDim(my.so, reduction = "pca")

# Create an expression matrix
expression_matrix <- my.so@assays$RNA@counts

# Get cell metadata
cell_metadata <- my.so@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}

# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(my.so@assays$RNA), row.names = rownames(my.so@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

# Seurat-derived CDS
my.cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)

# Transfer Seurat embeddings
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
reducedDim(my.cds, type = "PCA") <- my.so@reductions$pca@cell.embeddings 
my.cds@preprocess_aux$prop_var_expl <- my.so@reductions$pca@stdev
plot_pc_variance_explained(my.cds)

# Transfer Seurat UMAP embeddings
my.cds@int_colData@listData$reducedDims$UMAP <- my.so@reductions$umap@cell.embeddings
#    plot_cells(my.cds)

# Copy cluster info from Seurat
my.cds@clusters$UMAP_so$clusters <- my.so@meta.data$integrated_snn_res.0.3

my.cds <- cluster_cells(my.cds, reduction_method = "UMAP", resolution = 1e-3)

# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(my.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(my.cds@int_colData@listData$reducedDims$UMAP) <- NULL

DimPlot(my.so, reduction = "umap")
plot_cells(my.cds, group_label_size = 3.5)
plot_cells(my.cds, color_cells_by = "integrated_snn_res.0.3",
           show_trajectory_graph = T, group_label_size = 3.5)

my.cds <- learn_graph(my.cds, close_loop = F,learn_graph_control = list(minimal_branch_len=15, ncenter =500))

plot_cells(my.cds, color_cells_by = "integrated_snn_res.0.3",
           show_trajectory_graph = T, group_label_size = 3.5)


my.cds <- order_cells(my.cds)

###save it###
saveRDS(my.cds, file="../clustered_monocle.rds")

###### line plots to reduce drag of zeros 

library(RColorBrewer)
p<-brewer.pal(11,"Paired")


#make bin plots for figure 4
seur<-readRDS(file="../only_faps_integrated.rds")

my.cds<-readRDS(file="../clustered_monocle.rds")

#for each path
cds_sub_branch <- choose_graph_segments(my.cds, clear_cds = F)
#do below

pseudo_vals=pseudotime(cds_sub_branch)
rm(cds_sub_branch)
rm(my.cds)
rm(sub)

sub=subset(seur, cells=names(pseudo_vals))
identical(rownames(sub@meta.data),names(pseudo_vals))
sub$pseudo <- pseudo_vals
sub$Pseudo_bin <- as.numeric(cut_interval(sub$pseudo,20))
gene<-c("OSR1","GLI1","THY1","CD74","DLK1","MME","PIEZO2","SCX","DPP4","CD55","LOXL1","OSR2","VCAM1","SEMA3C","MSTN","TNMD","HGF","PTCH1"
        ,"LRIG1","ITGA8")
avg_exp=AverageExpression(sub, assay="RNA",features=gene,group.by = "Pseudo_bin")

#save them as you go

plot_val_fibro = melt(avg_exp)
plot_val_adipo = melt(avg_exp)
plot_val_teno = melt(avg_exp)

write.csv(plot_val_fibro, file=paste ("../plot_val_fibro.csv"))
write.csv(plot_val_adipo, file=paste ("../plot_val_adipo.csv"))
write.csv(plot_val_teno, file=paste ("../plot_val_teno.csv"))



##
plot_val_fibro_cd55 = subset(plot_val_fibro,Var1 %in% "CD55")
plot_val_adipo_cd55 = subset(plot_val_adipo,Var1 %in% "CD55")
plot_val_teno_cd55 = subset(plot_val_teno,Var1 %in% "CD55")

plot_val_fibro_cd55$group<-"fibro"
plot_val_adipo_cd55$group<-"adipo"
plot_val_teno_cd55$group<-"teno"

merged<-Reduce(function(x, y) merge(x, y, all=TRUE), 
               list(plot_val_fibro_cd55,
                    plot_val_adipo_cd55,
                    plot_val_teno_cd55))  

#write a loop
library(scales)

gene<-c("OSR1","GLI1","THY1","CD74","DLK1","MME","PIEZO2","SCX","DPP4","CD55")
gene<-c("LOXL1","OSR2","VCAM1","SEMA3C","MSTN","TNMD","HGF","PTCH1"
        ,"LRIG1","ITGA8")
gene<-c("MME","CD55")
for (i in 1:2){
  
  plot_val_fibro_gene<-subset(plot_val_fibro,Var1 %in% gene[i])
  plot_val_adipo_gene<-subset(plot_val_adipo,Var1 %in% gene[i])
  plot_val_teno_gene<-subset(plot_val_teno,Var1 %in% gene[i])
  
  
  p1=ggplot(plot_val_adipo_gene, aes(x=Var2,y=value,col=Var1,group=Var1))+
    geom_line(size=0.5,stat="identity",color="black")+
    geom_point(color="black")+
    theme_classic()+
    xlab("Pseudotime Bin")+
    ylab(paste0(gene[i],"expression"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          axis.text=element_text(size=20),axis.title=element_text(size=20))
  p2=ggplot(plot_val_fibro_gene, aes(x=Var2,y=value,col=Var1,group=Var1))+
    geom_line(size=0.5,stat="identity",color="black")+
    geom_point(color="black")+
    theme_classic()+
    xlab("Pseudotime Bin")+
    ylab("")+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          axis.text=element_text(size=20),axis.title=element_text(size=20))
  p3=ggplot(plot_val_teno_gene, aes(x=Var2,y=value,col=Var1,group=Var1))+
    geom_line(size=0.5,stat="identity",color="black")+
    geom_point(color="black")+
    theme_classic()+
    xlab("Pseudotime Bin")+
    ylab("")+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          axis.text=element_text(size=20),axis.title=element_text(size=20))
  
  ls<-c(max(plot_val_fibro_gene$value),
        max(plot_val_adipo_gene$value),
        max(plot_val_teno_gene$value))
  max(unlist(ls))
  
  pdf(paste0("../traj_log_lines_",gene[i],".pdf"),height=3,width=9)
 print(p1+scale_y_continuous(labels = label_number(accuracy = 1)
 )+coord_cartesian(ylim = c(0, max(unlist(ls)))) +
    p2+scale_y_continuous(labels = label_number(accuracy = 1)
    )+coord_cartesian(ylim = c(0, max(unlist(ls))))+
    p3+scale_y_continuous(labels = label_number(accuracy = 1)
    )+coord_cartesian(ylim = c(0, max(unlist(ls))))) 
  dev.off()
  
}
  
  ##

p1=ggplot(plot_val_fibro_cd55, aes(x=Var2,y=value,col=group,group=group))+
  geom_line(size=0.5,stat="identity",color="black")+
  geom_point(color="black")+
  theme_classic()+
  xlab("Pseudotime Bin")+
  ylab("CD55 expression")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p2=ggplot(plot_val_adipo_cd55, aes(x=Var2,y=value,col=group,group=group))+
  geom_line(size=0.5,stat="identity",color="black")+
  geom_point(color="black")+
  theme_classic()+
  xlab("Pseudotime Bin")+
  ylab("")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p3=ggplot(plot_val_teno_cd55, aes(x=Var2,y=value,col=group,group=group))+
  geom_line(size=0.5,stat="identity",color="black")+
  geom_point(color="black")+
  theme_classic()+
  xlab("Pseudotime Bin")+
  ylab("")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

pdf("../cd55_traj_log_lines.pdf",height=3,width=9)
p2+coord_cartesian(ylim = c(0, 20)) +
  p1+coord_cartesian(ylim = c(0, 20))+
  p3+coord_cartesian(ylim = c(0, 20)) 
dev.off()


#

###### make cell jitters for paths

gene<-c("OSR1","GLI1","THY1","CD74","DLK1","MME","PIEZO2","SCX","DPP4","CD55")
gene<-c("CD55")
#cols for pop 7 traj
fig1_col = c(p[1],p[2],p[3],p[5],p[6],p[7],p[8],p[9],p[10],p[11])
#cols for pop5  traj
fig1_col = c(p[1],p[2],p[3],p[5],p[6],p[7],p[8],p[9],p[11])
#cols for pop3  traj
fig1_col = c(p[1],p[2],p[3],p[4],p[5],p[7],p[8],p[9],p[11])

cds_sub_branch <- choose_graph_segments(my.cds, clear_cds = F)

pl=plot_genes_in_pseudotime(cds_sub_branch[gene,],
                            color_cells_by = "integrated_snn_res.0.3",panel_order = gene) 

dat=pl$data
dat$f_id <- factor(dat$f_id, levels =gene)

dat$integrated_snn_res.0.3<-factor(dat$integrated_snn_res.0.3, 
                                   levels=c("0","1",
                                            "2","3",
                                            "4","5",
                                            "6","7",
                                            "8","9",
                                            "10"))


#jitters only
pdf("../pop5_traj_log_jitter_1.5.pdf",height=1.5,width=8)
p2= ggplot(dat, aes(x=pseudotime,y=expectation))+
  rasterise(geom_jitter(aes(x=pseudotime,y=max(0),col=integrated_snn_res.0.3),size=0.1,width = 0.5,height = 0.1))+scale_color_manual(values=fig1_col)
p2+theme_void()+NoLegend()
dev.off()




