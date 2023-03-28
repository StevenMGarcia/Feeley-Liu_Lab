

import scvelo as scv
import pandas
from anndata import AnnData
import sys

#folder=sys.argv[1]
#print(folder)
#version=sys.argv[2]

#base='/myfiles/'

## import cellnames from seurat object
seur_names = scv.load("../cellnames/v1_cellnames.csv", index_col=0)

## import individual loom files in order to accurately match names with integrated Seurat object
adata = scv.read("../loomfiles/combined_v1_faps.loom")

adata.var_names_make_unique()
#clean names only for individual runs#
#scv.utils.clean_obs_names(adata)
adata=adata[seur_names["x"]]

scv.pp.filter_and_normalize(adata,min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

## dynamic modeling
##scv.tl.recover_dynamics(adata,var_names='all')
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

adata.write('/../loomfiles/v1_faps_veloc_dynamic.h5ad', compression='gzip')
