#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 17:52:53 2024

@author: maureen
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix


################################################################################################################################

# Import data
adata = sc.read_csv('/Users/maureen/Documents/projects/harris-lab/Isaac/hunter_dataset/GSM6280797_counts.csv').T
sc.pl.highest_expr_genes(adata, n_top=20)

################################################################################################################################

# PRE-PROCESSING

## annotate mitochondrial and ribosomal genes
adata.var['mt'] = adata.var.index.str.startswith('mt-')
adata.var['ribo'] = adata.var.index.str.startswith("Rps","Rpl")

## calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], jitter = 0.4, multi_panel = True) 
### "raw" data already appears pre-processed
### ribosomal reads are on scale for tissue type (https://kb.10xgenomics.com/hc/en-us/articles/218169723-What-fraction-of-reads-map-to-ribosomal-proteins-)

## standard filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

## filtering based on parameters in hunter paper (apparently already performed on raw data)
sc.pp.filter_cells(adata, min_counts = 1000)
sc.pp.filter_cells(adata, max_counts = 25000)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], jitter = 0.4, multi_panel = True) 

## slice adata object to perform actual filtering
adata = adata[adata.obs.total_counts > 1000]
adata = adata[adata.obs.total_counts < 25000]

adata = adata[adata.obs.pct_counts_mt < 20] 
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], jitter = 0.4, multi_panel = True) 


################################################################################################################################

# NORMALIZATION (LIBRARY CORRECTION)
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum = 1e4)
adata.layers['normalized'] = adata.X.copy()
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()
adata.raw = adata

# CLUSTERING
## identify highly variable genes (HVGs) that best describe the data
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable] ## filter out non HVGs
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo']) # helps remove sequencing artifacts

sc.pp.scale(adata, max_value = 10) ### Normalize each gene to unit variance of that gene
adata.layers['scaled'] = adata.X.copy()

sc.tl.pca(adata, svd_solver = 'arpack')  ### Run PCA for dimensionality reduction
sc.pl.pca_variance_ratio(adata, log = True, n_pcs = 50) #pick 30 given the elbow

## compute and embed neighborhood graph with the appropriate number of pcs
sc.pp.neighbors(adata, n_pcs = 30) ## calculate neighbors of cells based on top 30 PCs
sc.tl.umap(adata) # Embed and plot the UMAP
sc.pl.umap(adata)
sc.tl.leiden(adata, resolution = 0.7) #start with .5 --  which could not resolve CD4 vs CD8
sc.pl.umap(adata, color = ['leiden']) ## replot umap and color cells baed on leiden group


################################################################################################################################


# FIND CELL MARKERS AND IDENTIFY CLUSTERS

## Get markers
sc.tl.rank_genes_groups(adata, 'leiden', method = 'wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]

### Immune
sc.pl.umap(adata, color = ['leiden', 'Ptprc'], legend_loc = "on data") #CD45

#### Macrophages
sc.pl.umap(adata, color = ['leiden', 'Itgam'], legend_loc = "on data") #CD11b
sc.pl.umap(adata, color = ['leiden', 'Adgre1'], legend_loc = "on data") #PU.1
sc.pl.umap(adata, color = ['leiden', 'Hmox1'], legend_loc = "on data")#red pulp macrophages
sc.pl.umap(adata, color = ['leiden', 'Cd68'], legend_loc = "on data")


#### Monocytes
sc.pl.umap(adata, color = ['leiden', 'Itgam'], legend_loc = "on data") #CD11b
sc.pl.umap(adata, color = ['leiden', 'Ccr2'], legend_loc = "on data") 
sc.pl.umap(adata, color = ['leiden', 'Cmklr1'], legend_loc = "on data") 
sc.pl.umap(adata, color = ['leiden', 'Aif1'], legend_loc = "on data") 

#### Neutrophils
sc.pl.umap(adata, color = ['leiden', 'Itgam'], legend_loc = "on data") #CD11b
sc.pl.umap(adata, color = ['leiden', 'Ly6g'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Mmp9'], legend_loc = "on data") 
sc.pl.umap(adata, color = ['leiden', 'Csf3r'], legend_loc = "on data") 

#### Dendritic cells
sc.pl.umap(adata, color = ['leiden', 'Itgax'], legend_loc = "on data")  #DCs
sc.pl.umap(adata, color = ['leiden', 'Batf3'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Irf8'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Zbtb46'], legend_loc = "on data")  #DCs
sc.pl.umap(adata, color = ['leiden', 'Tlr9'], legend_loc = "on data")  #DCs


### T Cells / NK Cells
sc.pl.umap(adata, color = ['leiden', 'Cd3e'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Cd4'], legend_loc = "on data") 
sc.pl.umap(adata, color = ['leiden', 'Cd8a'], legend_loc = "on data") 
sc.pl.umap(adata, color = ['leiden', 'Tbx21'], legend_loc = "on data") 
sc.pl.umap(adata, color = ['leiden', 'Klre1'], legend_loc = "on data") #NK
sc.pl.umap(adata, color = ['leiden', 'Nkg7'], legend_loc = "on data") #NK

### B cells
sc.pl.umap(adata, color = ['leiden', 'Pxk'], legend_loc = "on data") #B cells
sc.pl.umap(adata, color = ['leiden', 'Cd19'], legend_loc = "on data") #B cells
sc.pl.umap(adata, color = ['leiden', 'Ighm'], legend_loc = "on data") #Marginal zone B cells
sc.pl.umap(adata, color = ['leiden', 'Mzb1'], legend_loc = "on data") # marignal zone B cell

### Erythroid lineage
sc.pl.umap(adata, color = ['leiden', 'Tfrc'], legend_loc = "on data") #erythroblasts
sc.pl.umap(adata, color = ['leiden', 'Gata1'], legend_loc = "on data") #erythroblasts
sc.pl.umap(adata, color = ['leiden', 'Tal1'], legend_loc = "on data") #erythroblasts
sc.pl.umap(adata, color = ['leiden', 'Hba-a1'], legend_loc = "on data") #erythroid-like
sc.pl.umap(adata, color = ['leiden', 'Hbb-bt'], legend_loc = "on data") #erythroid-like
sc.pl.umap(adata, color = ['leiden', 'Ermap'], legend_loc = "on data") #erythroblasts

### Endothelial cells
sc.pl.umap(adata, color = ['leiden', 'Pecam1'], legend_loc = "on data") 
sc.pl.umap(adata, color = ['leiden', 'Cdh5'], legend_loc = "on data") 
sc.pl.umap(adata, color = ['leiden', 'Ly6c1'], legend_loc = "on data") 


## Explore clusters 
markers[markers.group == '0'].head(20) # myeloid
markers[markers.group == '1'].head(20)
markers[markers.group == '2'].head(20)
markers[markers.group == '3'].head(20) # myeloid
markers[markers.group == '4'].head(20) #inflammatory
markers[markers.group == '5'].head(20) #myeloid, macrophage
markers[markers.group == '6'].head(20) #B cell
markers[markers.group == '7'].head(20) #myeloid
markers[markers.group == '8'].head(20) #T cell
markers[markers.group == '9'].head(20)
markers[markers.group == '10'].head(20) #NK cells / CD8 T cells
markers[markers.group == '11'].head(20) #possible monocytes
markers[markers.group == '12'].head(20)
markers[markers.group == '13'].head(20)
markers[markers.group == '14'].head(20) #YM1+ Macrophages
markers[markers.group == '15'].head(20) #Erythroid
markers[markers.group == '16'].head(20) #high ribosomal content
markers[markers.group == '17'].head(20) #myeloid
markers[markers.group == '18'].head(20) #Endothelial



## Create a dictionary to map cell labels
cell_type = {
"0":"B cells",
"1":"Erythroid-lineage",
"2":"Neutrophils / Granulocytes",
"3":"Ly6clo Monocytes",
"4":"Erythroid-lineage",
"5":"Marginal zone B cells / Plasmablasts",
"6":"CD4 T cells",
"7":"Neutrophils / Granulocytes",
"8":"Erythroid-lineage",
"9":"Macrophages",
"10":"CD8+ T cells / NK cells",
"11":"Ly6chi Monocytes",
"12":"Erythroid-lineage",
"13":"Neutrophils / Granulocytes",
"14":"NK cells / T cells",
"15":"Dendritic cells",
"16":"Erythroid-lineage",
"17":"Endothelial cells"
}

adata.obs['cell_type'] = adata.obs.leiden.map(cell_type)

sc.pl.umap(adata, color = 'cell_type')

################################################################################################################################

# DATA VISUALIZATION

## UMAP with cell type labels
adata.obs['cell type'] = adata.obs.leiden.map(cell_type)
sc.pl.umap(adata, color = ['cell type'])

## Set up dictionary for marker genes
marker_genes_dict = {
    'Hematopoietic': ['Ptprc'],
    'Mono / Macro': ['Adgre1', 'Hmox1', 'Cd68', 'Aif1', 'Cx3cr1'],
    'Neutrophils / Granulocytes': ['Ly6g', 'Mmp9', 'Csf3r'],
    'Dendritic cells': ['Itgax', 'Batf3', 'Zbtb46', 'H2-Oa'],
    'T cells': ['Cd3e', 'Cd4', 'Cd8a', 'Cd8b1'],
    'B cells': ['Cd19', 'Pxk', 'Ighm', 'Mzb1'],
    'Endothelial': ['Pecam1', 'Cdh5'],
    'Erythroid': ['Tfrc', 'Gata1', 'Hba-a1']
}



sc.pl.dotplot(adata, marker_genes_dict, groupby='cell_type', dendrogram = True)


################################################################################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/harris-lab/Isaac/hunter_dataset/write"

adata.X = csr_matrix(adata.X)

adata.write_h5ad(os.path.join(save_dir, 'clark-dataset-anotated.h5ad'))

################################################################################################################################






