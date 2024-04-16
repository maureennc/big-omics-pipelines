#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 13:36:34 2024

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

adata_concat = sc.read_h5ad(os.path.join(data_dir, 'concat-model_B.h5ad'))
adata_infected = sc.read_h5ad(os.path.join(data_dir, 'infected-only-model_B.h5ad'))

################################################################################################################################

## Import model_B
scvi_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/scvi'

model_B_dir = os.path.join(scvi_dir, 'model_B')
print(model_B_dir)
model_B = scvi.model.SCVI.load(model_B_dir, adata=adata_concat)
model_B

## Import model_C
model_C_dir = os.path.join(scvi_dir, 'model_C')
print(model_C_dir)
model_C = scvi.model.SCVI.load(model_C_dir, adata=adata_infected)
model_C

################################################################################################################################

# ANNOTATIONS

immune_cell_types = ['Macrophage', 'Microglia', 'CD4+ T cell', 'CD8+ T cell']
adata_infected.obs['immune_cell'] = adata_infected.obs['cell_type'].isin(immune_cell_types)

infiltrating = ['Macrophage', 'CD4+ T cell', 'CD8+ T cell']
adata_infected.obs['infiltrating'] = adata_infected.obs['cell_type'].isin(infiltrating)

t_cell_types = ['CD4+ T cell', 'CD8+ T cell']
adata_infected.obs['t_cell'] = adata_infected.obs['cell_type'].isin(t_cell_types)


adata_infected.obs['inflammatory_foci'] = adata_infected.obs['inflammatory_foci'].astype('category')

################################################################################################################################

# DE - BASE SCANPY

## Psuedobulk
sc.tl.rank_genes_groups(adata_infected, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(adata_infected, n_genes = 30, fontsize = 12)

## Infiltrating
infiltrating = adata_infected[(adata_infected.obs['infiltrating'] == True)].copy()

sc.tl.rank_genes_groups(infiltrating, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(infiltrating, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(infiltrating, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_infiltrating = markers[(markers.group == 'True')]
markers_infiltrating


## t_cell
t_cell = adata_infected[(adata_infected.obs['t_cell'] == True)].copy()

sc.tl.rank_genes_groups(t_cell, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(t_cell, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(t_cell, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_t_cell = markers[(markers.group == 'True')]
markers_t_cell



## Excitatory
excitatory = adata_infected[(adata_infected.obs['cell_type'] == 'Excitatory neuron')].copy()

sc.tl.rank_genes_groups(excitatory, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(excitatory, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(excitatory, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_excitatory = markers[(markers.group == 'True')]
markers_excitatory

## Inhibitory
inhibitory = adata_infected[(adata_infected.obs['cell_type'] == 'Inhibitory neuron')].copy()

sc.tl.rank_genes_groups(inhibitory, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(inhibitory, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(inhibitory, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_inhibitory = markers[(markers.group == 'True')]
markers_inhibitory


## Vascular
vascular = adata_infected[(adata_infected.obs['cell_type'] == 'Vascular')].copy()

sc.tl.rank_genes_groups(vascular, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(vascular, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(vascular, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_vascular = markers[(markers.group == 'True')]
markers_vascular


## Choroid plexus
cp = adata_infected[(adata_infected.obs['cell_type'] == 'Choroid plexus')].copy()

sc.tl.rank_genes_groups(cp, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(cp, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(vascular, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_cp = markers[(markers.group == 'True')]
markers_cp


## Astrocyte
astrocyte = adata_infected[(adata_infected.obs['cell_type'] == 'Astrocyte')].copy()

sc.tl.rank_genes_groups(astrocyte, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(astrocyte, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(astrocyte, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_astrocyte = markers[(markers.group == 'True')]
markers_astrocyte

## Microglia
microglia = adata_infected[(adata_infected.obs['cell_type'] == 'Microglia')].copy()

sc.tl.rank_genes_groups(microglia, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(microglia, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(microglia, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_microglia = markers[(markers.group == 'True')]
markers_microglia


## Macrophage
macrophage = adata_infected[(adata_infected.obs['cell_type'] == 'Macrophage')].copy()

sc.tl.rank_genes_groups(macrophage, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(macrophage, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(macrophage, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_mac = markers[(markers.group == 'True')]
markers_mac


## CD4+
cd4 = adata_infected[(adata_infected.obs['cell_type'] == 'CD4+ T cell')].copy()


sc.tl.rank_genes_groups(cd4, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(cd4, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(cd4, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_cd4 = markers[(markers.group == 'True')]
markers_cd4


## cd8+
cd8 = adata_infected[(adata_infected.obs['cell_type'] == 'CD8+ T cell')].copy()

sc.tl.rank_genes_groups(cd8, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(cd8, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(cd8, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_cd8 = markers[(markers.group == 'True')]
markers_cd8

## Oligodendrocyte
oligo = adata_infected[(adata_infected.obs['cell_type'] == 'Oligodendrocyte')].copy()

sc.tl.rank_genes_groups(oligo, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(oligo, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(oligo, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_oligo = markers[(markers.group == 'True')]
markers_oligo

## OPC
opc = adata_infected[(adata_infected.obs['cell_type'] == 'OPC')].copy()

sc.tl.rank_genes_groups(opc, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(opc, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(opc, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_opc = markers[(markers.group == 'True')]
markers_opc


################################################################################################################################

## immune DE for VOLCANO PLOT
immune = adata_infected[(adata_infected.obs['immune_cell'] == True)].copy()

sc.tl.rank_genes_groups(immune, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(immune, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(immune, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_immune = markers[(markers.group == 'True')]
markers_immune

save_path = "/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/spreadsheets"
markers_immune.to_csv(os.path.join(save_path, 'cluster-vs-rest.csv'), index = False)



######

# HOTSPOT CD4 vs CD8
adata_infected = sc.read_h5ad(os.path.join(data_dir, 'infected-only-model_B.h5ad'))

tcells = adata_infected[adata_infected.obs['cell_type'].isin(['CD4+ T cell', 'CD8+ T cell'])].copy()
print(tcells.obs.cell_type.value_counts())
      
tcells.obs['inflammatory_foci'] = tcells.obs['inflammatory_foci'].astype('category')
tcells.obs['inflammatory_foci'] = tcells.obs['inflammatory_foci'].astype('str')

print(tcells.obs.inflammatory_foci.value_counts())
tcells = tcells[tcells.obs['inflammatory_foci'] == 'True'].copy()
tcells

sc.tl.rank_genes_groups(tcells, groupby='cell_type', method='wilcoxon')
sc.pl.rank_genes_groups(tcells, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(tcells, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_tcells = markers[(markers.group == 'CD8+ T cell')]
markers_tcells


## T cells in clusters vs rest
adata_infected = sc.read_h5ad(os.path.join(data_dir, 'infected-only-model_B.h5ad'))

tcells = adata_infected[adata_infected.obs['cell_type'].isin(['CD4+ T cell', 'CD8+ T cell'])].copy()
print(tcells.obs.cell_type.value_counts())
      
tcells.obs['inflammatory_foci'] = tcells.obs['inflammatory_foci'].astype('category')
tcells.obs['inflammatory_foci'] = tcells.obs['inflammatory_foci'].astype('str')


sc.tl.rank_genes_groups(tcells, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(tcells, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(tcells, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_tcells = markers[(markers.group == 'True')]
markers_tcells


##########################

# FOR POSTER


# Microlia clusters vs rest
## Microglia


microglia = adata_infected[(adata_infected.obs['cell_type'] == 'Microglia')].copy()
microglia.obs['inflammatory_foci'] = microglia.obs['inflammatory_foci'].astype('category')


sc.tl.rank_genes_groups(microglia, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(microglia, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(microglia, None)
markers = markers[(markers.pvals_adj < 0.05)]
markers_microglia = markers[(markers.group == 'True')]
markers_microglia

save_path = "/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/spreadsheets"
markers_microglia.to_csv(os.path.join(save_path, 'microglia-cluster-vs-rest.csv'), index = False)



# INFILTRATING CLUSTER VS REST
adata_infected = sc.read_h5ad(os.path.join(data_dir, 'infected-only-model_B.h5ad'))

infiltrating = adata_infected[adata_infected.obs['cell_type'].isin(['CD4+ T cell', 'CD8+ T cell', 'Macrophage'])].copy()
print(infiltrating.obs.cell_type.value_counts())
      
infiltrating.obs['inflammatory_foci'] = infiltrating.obs['inflammatory_foci'].astype('category')
infiltrating.obs['inflammatory_foci'] = infiltrating.obs['inflammatory_foci'].astype('str')

print(infiltrating.obs.inflammatory_foci.value_counts())



sc.tl.rank_genes_groups(infiltrating, groupby='inflammatory_foci', method='wilcoxon')
sc.pl.rank_genes_groups(infiltrating, n_genes = 30, fontsize = 12)

markers = sc.get.rank_genes_groups_df(infiltrating, None)
#markers = markers[(markers.pvals_adj < 0.05)]
markers_infiltrating = markers[(markers.group == 'True')]
markers_infiltrating


save_path = "/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/spreadsheets"
markers_infiltrating.to_csv(os.path.join(save_path, 'infiltrating-cluster-vs-rest.csv'), index = False)


################################################################################################################################
# DE SUMMARY


df_list = [markers_excitatory, markers_inhibitory, markers_vascular, markers_cp, markers_astrocyte, markers_microglia, markers_mac, markers_cd4, markers_cd8, markers_oligo, markers_opc] 
cell_types = ['Excitatory neuron', 'Inhibitory neuron', 'Vascular', 'Choroid plexus', 'Astrocyte', 'Microglia', 'Macrophage', 'CD4+ T cell', 'CD8+ T cell', 'Oligodendrocyte', 'OPC']

de_counts = {}

for cell_type, df in zip(cell_types, df_list):
    count_significant = (df['pvals_adj'] < 0.05).sum()  # or df['is_de_fdr_0.05'].sum()
    de_counts[cell_type] = count_significant

de_counts_df = pd.DataFrame(list(de_counts.items()), columns=['CellType', 'SignificantDECount'])
print(de_counts_df)

### Bar plot
plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
barplot = sns.barplot(x='CellType', y='SignificantDECount', data=de_counts_df, palette='viridis')

plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
plt.xlabel('Cell Type')  # X-axis label
plt.ylabel('Number of Significantly DE Genes')  # Y-axis label
plt.title('Significantly DE Genes')  # Plot title

# Annotating bars with the count of significantly DE genes
for p in barplot.patches:
    barplot.annotate(format(p.get_height(), '.0f'), 
                     (p.get_x() + p.get_width() / 2., p.get_height()), 
                     ha='center', va='center', 
                     xytext=(0, 9), 
                     textcoords='offset points')

plt.tight_layout()
plt.show()

################################################################################################################################


# DATA VISUALIZATION - CLUSTER VS. REST


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
mg = ['Cxcl9', 'Cxcl10', 'C3', 'Cybb', 'Tnf', 'Vcam1', 'Tnfaip2', 'Il12rb1', 'Clec7a', 'Itgal', 'Icam1', 'Ccl2', 'Itgax', 'Socs1', 'Axl', 'Il1a', 'Il1b', 'P2ry12', 'Tmem119', 'Trem2', 'Sall1', 'Adora3', 'Csf1r', 'Cx3cr1', 'Syk', 'Stat1', 'Nlrp1a', 'Nlrp1b', 'S100a9']


sc.pl.matrixplot(adata_infected, cell_death_genes, groupby = 'inflammatory_foci', dendrogram = True, use_raw = True, cmap = 'bwr', standard_scale = 'var') #standard_scale ='var'
sc.pl.matrixplot(adata_infected, cytokine, groupby = 'inflammatory_foci', dendrogram = True, use_raw = True,  cmap = 'bwr', standard_scale = 'var')
sc.pl.matrixplot(adata_infected, chemokine, groupby = 'inflammatory_foci', dendrogram = True, use_raw = True,  cmap = 'bwr', standard_scale = 'var')
sc.pl.matrixplot(adata_infected, gf, groupby = 'inflammatory_foci', dendrogram = True, use_raw = True,  cmap = 'bwr', standard_scale = 'var')
sc.pl.matrixplot(adata_infected, cell_activation, groupby = 'inflammatory_foci', dendrogram = True, use_raw = True,  cmap = 'bwr', standard_scale = 'var')
sc.pl.matrixplot(adata_infected, dam_activation, groupby = 'inflammatory_foci', dendrogram = True, use_raw = True,  cmap = 'bwr', standard_scale = 'var')

sc.pl.matrixplot(adata_infected, mg, groupby = 'inflammatory_foci', dendrogram = True, use_raw = True,  cmap = 'bwr', standard_scale = 'var')

genes = ['Ifngr1',  'Stat1', 'Socs1', 'Tap1', 'Icam1', 'Vcam1', 'Ccl2', 'Cxcl9', 'Cxcl10', 'Rel', 'Rela', 'Relb', 'Nfkb1', 'Zbp1', 'Ripk1', 'Ripk3', 'Casp3', 'Casp8']
sc.pl.matrixplot(adata_infected, genes, groupby = 'inflammatory_foci', dendrogram = True, use_raw = True,  cmap = 'viridis', standard_scale = 'var')


################################################################################################################################


# CLUSTER VS. REST - scVI

## Immune
idx1 = (adata_infected.obs['immune_cell'] == True) & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['immune_cell'] == True) & (adata_infected.obs['inflammatory_foci'] == False)

pseudobulk_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

pseudobulk_de

## Infiltrating
idx1 = (adata_infected.obs['infiltrating'] == True) & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['infiltrating'] == True) & (adata_infected.obs['inflammatory_foci'] == False)

infiltrating_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

infiltrating_de

## T cell
idx1 = (adata_infected.obs['t_cell'] == True) & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['t_cell'] == True) & (adata_infected.obs['inflammatory_foci'] == False)

tcell_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

tcell_de

## Astrocyte
idx1 = (adata_infected.obs['cell_type'] == 'Astrocyte') & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['cell_type'] == 'Astrocyte') & (adata_infected.obs['inflammatory_foci'] == False)

astrocyte_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

astrocyte_de


## Microglia
idx1 = (adata_infected.obs['cell_type'] == 'Microglia') & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['cell_type'] == 'Microglia') & (adata_infected.obs['inflammatory_foci'] == False)

microglia_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

microglia_de

## Macrophage
idx1 = (adata_infected.obs['cell_type'] == 'Macrophage') & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['cell_type'] == 'Macrophage') & (adata_infected.obs['inflammatory_foci'] == False)

macrophage_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

macrophage_de


## Vascular
idx1 = (adata_infected.obs['cell_type'] == 'Vascular') & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['cell_type'] == 'Vascular') & (adata_infected.obs['inflammatory_foci'] == False)

vascular_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

vascular_de


## CD4+ T cell
idx1 = (adata_infected.obs['cell_type'] == 'CD4+ T cell') & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['cell_type'] == 'CD4+ T cell') & (adata_infected.obs['inflammatory_foci'] == False)

cd4_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

cd4_de


## CD8+ T cell
idx1 = (adata_infected.obs['cell_type'] == 'CD8+ T cell') & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['cell_type'] == 'CD8+ T cell') & (adata_infected.obs['inflammatory_foci'] == False)

cd8_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

cd8_de

## Excitatory neuron
idx1 = (adata_infected.obs['cell_type'] == 'Excitatory neuron') & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['cell_type'] == 'Excitatory neuron') & (adata_infected.obs['inflammatory_foci'] == False)

excitatory_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

excitatory_de

## Inhibitory neuron
idx1 = (adata_infected.obs['cell_type'] == 'Inhibitory neuron') & (adata_infected.obs['inflammatory_foci'] == True)
idx2 = (adata_infected.obs['cell_type'] == 'Inhibitory neuron') & (adata_infected.obs['inflammatory_foci'] == False)

Inhibitory_de = model_C.differential_expression(
    idx1 = idx1,
    idx2 = idx2)

Inhibitory_de


################################################################################################################################



# RUN DE WITH BASE SCANPY

subset = adata_infected[(adata_infected.obs['condition'] == 'treated') & (adata_infected.obs['cell_type'] == 'neuron')].copy()

sc.tl.rank_genes_groups(adata_temp, groupby='leiden_scVI', method='wilcoxon')


################################################################################################################################


################################################################################################################################
################################################################################################################################

# NOTES

