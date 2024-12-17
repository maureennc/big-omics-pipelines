#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:44:05 2024

@author: maureen
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

import random
import torch
import scvi

print(sns.__version__)
print(pd.__version__)
print(np.__version__)
print(sc.__version__)
print(scvi.__version__)

################################################################################################################################

# SETTINGS

## Random seed
random.seed(0)
torch.manual_seed(0)
np.random.seed(0)
scvi.settings.seed = 0

## Matplotlib
%matplotlib qt5
plt.rcParams['font.family'] = 'Arial'

################################################################################################################################

# IMPORT DATA

## H5ad
data_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/h5ad/processed-data'

adata = sc.read_h5ad(os.path.join(data_dir, 'naive-model_A-roi-annotated.h5ad'))

################################################################################################################################

# GENE LISTS

cell_death_genes = ['Fasl', 'Ager', 'Casp8', 'Casp3', 'Casp7', 'Casp9', 'Tnfrsf1a', 'Tnfrsf1b', 'Ripk1', 'Ripk3', 'Fas', 
                    'Tnfrsf10b', 'Tnfrsf25', 'Tnfrsf21', 'Fadd', 'Tradd', 'Apaf1', 'Cflar', 'Birc2', 'Birc3', 'Mlkl', 'Cyld', 
                    'Traf2', 'Chuk', 'Ikbkg', 'Nfkbia', 'Tlr4', 'Ticam1', 'Zbp1', 'Bax', 'Bak1', 'Bad', 'Bid', 'Bcl2l11', 
                    'Pmaip1', 'Bbc3', 'Bmf', 'Blk', 'Casp1', 'Casp4', 'Nlrp3', 'Nlrc4', 'Nlrc3', 'Nlrp1a', 'Nlrp1b', 'Nlrp1c-ps',
                    'Nlrp12', 'Naip6', 'Naip2', 'Aim2', 'Gsdmd', 'Gsdme', 'Il1a', 'Dpp8', 'Dpp9', 'Nek7', 'Pycard', 
                    'S100a9', 'Tnf']

gas6_genes = ['Gas6', 'Ager', 'Ntrk2', 'Tgfbr2', 'Il6st', 'Vcam1', 'Cdh5', 'Stat2', 'Cflar', 'Bad', 'Tap1', 'C3',
              'Tgfb2', 'Aldh1l1', 'Tnfrsf1a', 'Vegfa', 'Vegfb', 'Vegfc', 'Stat1', 'Il17rc', 'Zbp1', 'Mpeg1', 'Il10rb', 'Rorc', 'Csf1',
              'Nfkb1', 'Traf2', 'Ifngr1', 'Stat1', 'Ifnar2', 'Dtna', 'Nes', 'Cxcl10', 'Fas', 'Bbc3', 'Ikbkg', 'Rela', 'Nfkbia', 
              'Zbp1', 'Casp8']

cytokine = ['Il1a', 'Il1b', 'Il1r1', 'Il1r2', 'Il2', 'Il2ra', 'Il2rb', 'Il2rg', 'Il4', 'Il4ra', 'Il6', 'Il6ra', 'Il7', 'Il7r', 'Il9', 'Il9r', 'Il10', 'Il10ra', 'Il10rb', 'Il12a', 'Il12b', 'Il12rb1', 'Il12rb2', 'Il13', 'Il13ra1', 'Il15', 'Il15ra', 'Il17a', 'Il17b', 'Il17c', 'Il17d', 'Il17f', 'Il17ra', 'Il17rb', 'Il17rc', 'Il18', 'Il21', 'Il21r', 'Il23a', 'Il23r', 'Il27', 'Il27ra', 'Il6st', 'Il33', 'Il1rl1', 'Il1rap', 'Il34', 'Csf1r', 'Ebi3', 'Il1rl2', 'Ifnab', 'Ifnb1', 'Ifnar1', 'Ifnar2', 'Ifng', 'Ifngr1', 'Ifngr2', 'Ifnlr1', 'Myd88', 'Socs1']
chemokine = ['Ccl2', 'Ccl19', 'Ccr7', 'Ccr5', 'Cxcl1', 'Cxcl2', 'Cxcr2', 'Cxcl9', 'Cxcl10', 'Cxcr3', 'Cx3cl1', 'Cxcr6', 'Cxcl16', 'Xcr1']
gf = ['Csf1', 'Csf2ra', 'Csf2rb', 'Csf2rb2', 'Csf3', 'Csf3r', 'Vegfa', 'Vegfb', 'Vegfc', 'Vegfd', 'Flt4', 'Tgfb1', 'Tgfb2', 'Tgfb3', 'Tgfbr1', 'Bdnf', 'Ntrk2', 'Ngf', 'Ntrk1']
cell_activation = ['Itgal', 'Cd44', 'Itga1', 'Itga4', 'Vcam1', 'Icam1', 'Prf1', 'Tap1', 'Mmp9', 'Cd69', 'Cd80', 'Cd86', 'Klrg1', 'Fos', 'Fosb', 'Jun', 'Nr4a1', 'Myc', 'Egr1', 'Arc', 'Nfatc1', 'Nfatc2', 'Nfatc3', 'Nfatc4', 'Nfat5', 'Nfkb1', 'Nfkb2', 'Rel', 'Rela', 'Relb', 'Stat1', 'Stat2', 'Stat3', 'Stat4', 'Stat5a', 'Stat5b', 'Stat6']
cell_death = ['Fasl', 'Ager', 'Casp8', 'Casp3', 'Casp7', 'Casp9', 'Tnfrsf1a', 'Tnfrsf1b', 'Ripk1', 'Ripk3', 'Fas', 'Tnfrsf10b', 'Tnfrsf25', 'Tnfrsf21', 'Fadd', 'Tradd', 'Apaf1', 'Cflar', 'Birc2', 'Birc3', 'Mlkl', 'Cyld', 'Traf2', 'Chuk', 'Ikbkg', 'Nfkbia', 'Tlr4', 'Ticam1', 'Zbp1', 'Bax', 'Bak1', 'Bad', 'Bid', 'Bcl2l11', 'Pmaip1', 'Bbc3', 'Bmf', 'Blk', 'Hrk', 'Casp1', 'Casp4', 'Nlrp3', 'Nlrc4', 'Nlrc3', 'Nlrp1a', 'Nlrp1b', 'Nlrp12', 'Naip6', 'Naip2', 'Nlrp1c-ps', 'Aim2', 'Gsdma', 'Gsdmd', 'Gsdme', 'Dpp8', 'Dpp9', 'Nek7', 'Pycard', 'S100a9', 'Tnf', 'Il1a']
dam_activation = ['Axl', 'Msr1', 'Itgax', 'Spp1', 'Arg1', 'Siglec1', 'Gpnmb', 'Ccrl2', 'Lpl', 'Trem2', 'Clec7a', 'Syk', 'Mpeg1']

genes = ['Ifngr1',  'Stat1', 'Socs1', 'Tap1', 'Icam1', 'Vcam1', 'Ccl2', 'Cxcl9', 'Cxcl10', 'Rel', 'Rela', 'Relb', 'Nfkb1', 'Zbp1', 'Ripk1', 'Ripk3', 'Casp3', 'Casp8']

genes_df = pd.DataFrame(adata.var)
################################################################################################################################

# NAIVE DE - TESTING

sc.tl.rank_genes_groups(adata, groupby='region_of_interest', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers

sc.tl.dendrogram(adata, groupby = 'region_of_interest')

sc.pl.matrixplot(adata, cell_death, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(adata, gas6_genes, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(adata, cytokine, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(adata, chemokine, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(adata, gf, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(adata, cell_activation, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(adata, dam_activation, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)

# HERE

sc.pl.matrixplot(adata, genes, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True, cmap = 'inferno')

################################################################################################################################

# FIX REGION ANNOTATIONS (CAPITALIZE)
annotations_to_capitalize = ['cortex', 'hippocampus', 'thalamus', 'hypothalamus', 'striatum', 'ventricle']

def capitalize_annotation(annotation):
    if annotation in annotations_to_capitalize:
        return annotation.capitalize()
    else:
        return annotation
adata.obs['region_of_interest'] = adata.obs['region_of_interest'].apply(capitalize_annotation)

adata.obs['region_of_interest'].value_counts()

################################################################################################################################

# NAIVE DE SUMMARY

adata.obs['region_of_interest'] = adata.obs['region_of_interest'].astype('category')

sc.tl.rank_genes_groups(adata, groupby='region_of_interest', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers



region = ['Cortex', 'Hippocampus', 'Thalamus', 'Hypothalamus', 'Striatum', 'Ventricle']

de_counts = {}

for cell_type in region:
    # Filter the DataFrame for the current region
    filtered_markers = markers[markers['group'] == cell_type]
    # Count the number of significant DE genes for the filtered DataFrame
    count_significant = (filtered_markers['pvals_adj'] < 0.05).sum()
    de_counts[cell_type] = count_significant

# Convert the summary data into a DataFrame
de_counts_df = pd.DataFrame(list(de_counts.items()), columns=['Region', 'SignificantDECount'])

print(de_counts_df)


## Seaborn customizaiton

plt.figure(figsize=(5, 3), dpi=500)  
barplot = sns.barplot(x='Region', y='SignificantDECount', data=de_counts_df, palette='inferno')

plt.xticks(rotation=45, ha='right') 
plt.xlabel('')  # X-axis label
plt.ylabel('# of Significant\nGenes ')  
plt.title('')  # Plot title

# Annotating bars with the count of significantly DE genes
for p in barplot.patches:
    barplot.annotate(format(p.get_height(), '.0f'), 
                     (p.get_x() + p.get_width() / 2., p.get_height()), 
                     ha='center', va='center', 
                     xytext=(0, 6), 
                     fontsize=12, 
                     textcoords='offset points')
max_height = max(p.get_height() for p in barplot.patches)
plt.ylim(0, max_height + (max_height * 0.2))  # Add 20% headroom above the highest bar
plt.tight_layout()
plt.grid(False)
plt.show()


### Bar plot viz old
plt.figure(figsize=(5, 3), dpi=500)  
barplot = sns.barplot(x='Region', y='SignificantDECount', data=de_counts_df, palette='viridis')

plt.xticks(rotation=45, ha='right')  
plt.xlabel('')
plt.ylabel('# of Significant \nGenes ') 
plt.title('') 

# Annotating bars with the count of significantly DE genes
for p in barplot.patches:
    barplot.annotate(format(p.get_height(), '.0f'), 
                     (p.get_x() + p.get_width() / 2., p.get_height()), 
                     ha='center', va='center', 
                     xytext=(0, 6),
                     fontsize = 12,
                     textcoords='offset points')
max_height = max(p.get_height() for p in barplot.patches)
plt.ylim(0, max_height + (max_height * 0.2))  # Adding 20% headroom above the highest bar
plt.tight_layout()
plt.grid(False)
plt.show()




################################################################################################################################

# OPTIMIZED DE VIZ

region = ['Cortex', 'Hippocampus', 'Thalamus', 'Hypothalamus', 'Striatum', 'Ventricle']

de_counts = {}

for cell_type in region:
    filtered_markers = markers[markers['group'] == cell_type]
    count_significant = (filtered_markers['pvals_adj'] < 0.05).sum()
    de_counts[cell_type] = count_significant

de_counts_df = pd.DataFrame(list(de_counts.items()), columns=['Region', 'SignificantDECount'])
de_counts_df = de_counts_df.sort_values(by='SignificantDECount')

plt.figure(figsize=(5, 3), dpi=500) 
barplot = sns.barplot(x='Region', y='SignificantDECount', data=de_counts_df, palette='inferno')

plt.xticks(rotation=45, ha='right')  
plt.xlabel('') 
plt.ylabel('# of Significant\nGenes ') # Text wrapping
plt.title('')

# Annotating bars with the count of significantly DE genes
for p in barplot.patches:
    barplot.annotate(format(p.get_height(), '.0f'), 
                     (p.get_x() + p.get_width() / 2., p.get_height()), 
                     ha='center', va='center', 
                     xytext=(0, 6), 
                     fontsize=12, 
                     textcoords='offset points')

max_height = max(p.get_height() for p in barplot.patches)
plt.ylim(0, max_height + (max_height * 0.2))
plt.tight_layout()
plt.grid(False)
plt.show()




################################################################################################################################

# MICROGLIA DE

microglia = adata[(adata.obs['cell_type'] == 'Microglia')].copy()

sc.pl.matrixplot(microglia, cell_death, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(microglia, cytokine, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(microglia, chemokine, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(microglia, gf, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(microglia, cell_activation, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)


sc.tl.rank_genes_groups(microglia, groupby='region_of_interest', method='wilcoxon')
sc.pl.rank_genes_groups(microglia, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers

################################################################################################################################

# astrocyte DE

astrocyte = adata[(adata.obs['cell_type'] == 'Astrocyte')].copy()

sc.pl.matrixplot(astrocyte, cell_death, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(astrocyte, cytokine, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(astrocyte, chemokine, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(astrocyte, gf, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(astrocyte, cell_activation, groupby = 'region_of_interest', standard_scale = 'var', use_raw = True, dendrogram = True)


sc.tl.rank_genes_groups(astrocyte, groupby='region_of_interest', method='wilcoxon')
sc.pl.rank_genes_groups(astrocyte, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers
