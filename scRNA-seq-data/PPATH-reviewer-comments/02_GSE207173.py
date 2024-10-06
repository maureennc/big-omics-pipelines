#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 11:28:23 2024

@author: maureen
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
sc.set_figure_params(scanpy = True, dpi = 500, dpi_save = 1000)

import seaborn as sns

save_dir = "/Users/maureen/Documents/projects/harris-lab/Isaac/hunter_dataset/figures"
os.chdir(save_dir)
################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/harris-lab/Isaac/hunter_dataset/write"

adata = sc.read_h5ad(os.path.join(data_dir, 'clark-dataset-anotated.h5ad'))

adata.obs['cell_type'] = adata.obs['cell_type'].replace('CD4 T cells', 'CD4+ T cells')


################################################################################################################################

# ANNOTATE CELLS WITH NON-ZERO CASP1 TRANSCRIPT COUNT

casp1_counts = adata[:, 'Casp1'].layers['counts']
non_zero_casp1_mask = casp1_counts.X if hasattr(casp1_counts, 'X') else casp1_counts > 0
adata.obs['casp1_positive'] = non_zero_casp1_mask
print(adata.obs['casp1_positive'].value_counts())


adata.obs['casp1_expression'] = adata.obs['casp1_positive'].map({True: 'casp1-positive', False: 'casp1-negative'})
print(adata.obs['casp1_expression'].value_counts())

################################################################################################################################

# ANNOTATE CELLS WITH NON-ZERO CX3CR1 TRANSCRIPT COUNT

cx3cr1_counts = adata[:, 'Cx3cr1'].layers['counts']
non_zero_cx3cr1_mask = cx3cr1_counts.X if hasattr(cx3cr1_counts, 'X') else cx3cr1_counts > 0
adata.obs['Cx3cr1_positive'] = non_zero_cx3cr1_mask
print(adata.obs['Cx3cr1_positive'].value_counts())


adata.obs['Cx3cr1_expression'] = adata.obs['Cx3cr1_positive'].map({True: 'Cx3cr1-positive', False: 'Cx3cr1-negative'})
print(adata.obs['Cx3cr1_expression'].value_counts())


################################################################################################################################

# FIND INTERSECTION OF CX3CR1 AND CASP1 EXPRESSION

## Define a mask for Casp1+ Cx3cr1+ cells
double_positive_mask = (adata.obs['casp1_positive'] == True) & (adata.obs['Cx3cr1_positive'] == True)

## Define a mask for Casp1- Cx3cr1+ cells
single_positive_mask = (adata.obs['casp1_positive'] == False) & (adata.obs['Cx3cr1_positive'] == True)

# Annotate these groups in the dataframe
adata.obs['comparison_group'] = 'None'  # Initialize the column with 'None'
adata.obs.loc[double_positive_mask, 'comparison_group'] = 'Casp1+Cx3cr1+'
adata.obs.loc[single_positive_mask, 'comparison_group'] = 'Casp1-Cx3cr1+'


################################################################################################################################

# DIFFERENTIAL EXPRESSION (IMMUNE CELLS)

## NEED TO CHECK THIS SSECTION BELOW; MASKING INCORRECT
mdata = adata[adata.obs['cell_type'] != 'Erythroid-lineage'].copy()
mdata = mdata[mdata.obs['cell_type'] != 'Endothelial cells'].copy()
mdata = mdata[mdata.obs['comparison_group'] != 'None'].copy()

print(mdata.obs['comparison_group'].value_counts())


sc.tl.rank_genes_groups(mdata, groupby='comparison_group', method='wilcoxon', use_raw = False, layer = 'log1p')
sc.pl.rank_genes_groups(mdata, n_genes = 30, fontsize = 8)

de_results = sc.get.rank_genes_groups_df(mdata, None)
de_results = de_results[(de_results.group == 'Casp1+Cx3cr1+')].copy()

save_path = "/Users/maureen/Documents/projects/harris-lab/Isaac/hunter_dataset/write"

#de_results.to_csv(os.path.join(save_path, 'casp1-pos-vs-neg-DE-full.csv'), index = False)

de_results = de_results[(de_results.pvals_adj < 0.05)]

#de_results.to_csv(os.path.join(save_path, 'casp1-pos-vs-neg-DE-significant.csv'), index = False)

de_results_upregulated = de_results[(de_results.logfoldchanges > 0)]
de_results_downregulated = de_results[(de_results.logfoldchanges < 0)]

#de_results_upregulated.to_csv(os.path.join(save_path, 'casp1-pos-vs-neg-DE-up.csv'), index = False)
#de_results_downregulated.to_csv(os.path.join(save_path, 'casp1-pos-vs-neg-DE-down.csv'), index = False)

################################################################################################################################

# BAR PLOT

## Casp1+
casp1_positive_data = mdata[mdata.obs['comparison_group'] == 'Casp1+Cx3cr1+']

cell_type_counts = casp1_positive_data.obs['cell_type'].value_counts()
cell_type_frequencies = (cell_type_counts / cell_type_counts.sum()) * 100

cell_type_frequencies_df = cell_type_frequencies.reset_index()
cell_type_frequencies_df.columns = ['Cell Type', 'Frequency (%)']

order_list = ['CD8+ T cells / NK cells', 'CD4+ T cells', 'Macrophages', 'Ly6clo Monocytes', 'Ly6chi Monocytes', 'Dendritic cells', 'Marginal zone B cells / Plasmablasts', 'Neutrophils / Granulocytes', 'NK cells / T cells']

plt.figure(figsize=(6, 2), dpi = 500)
sns_barplot = sns.barplot(x='Cell Type', y='Frequency (%)', data=cell_type_frequencies_df, palette="rocket", order = order_list)
sns_barplot.set_xticklabels(sns_barplot.get_xticklabels(), rotation=45, ha='right')
sns_barplot.set_title('Distribution of Casp1-positive Immune Cells')
sns_barplot.set_xlabel('')
sns_barplot.set_ylabel('Frequency (%)')

plt.tight_layout()
plt.grid(False)
plt.show()

## Casp1-

casp1_positive_data = mdata[mdata.obs['comparison_group'] == 'Casp1-Cx3cr1+']

cell_type_counts = casp1_positive_data.obs['cell_type'].value_counts()
cell_type_frequencies = (cell_type_counts / cell_type_counts.sum()) * 100

cell_type_frequencies_df = cell_type_frequencies.reset_index()
cell_type_frequencies_df.columns = ['Cell Type', 'Frequency (%)']

order_list = ['CD8+ T cells / NK cells', 'CD4+ T cells', 'Macrophages', 'Ly6clo Monocytes', 'Ly6chi Monocytes', 'Dendritic cells', 'Marginal zone B cells / Plasmablasts', 'Neutrophils / Granulocytes', 'NK cells / T cells']

plt.figure(figsize=(6, 2), dpi = 500)
sns_barplot = sns.barplot(x='Cell Type', y='Frequency (%)', data=cell_type_frequencies_df, palette="rocket", order = order_list)
sns_barplot.set_xticklabels(sns_barplot.get_xticklabels(), rotation=45, ha='right')
sns_barplot.set_title('Distribution of Casp1-negative Immune Cells')
sns_barplot.set_xlabel('')
sns_barplot.set_ylabel('Frequency (%)')

plt.tight_layout()
plt.grid(False)
plt.show()

################################################################################################################################

# GROUPED BAR PLOT

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Function to calculate frequency percentages
def get_frequency_percentages(data, group_name):
    # Calculate the frequency of each cell type within the group
    cell_type_counts = data.obs['cell_type'].value_counts()
    cell_type_frequencies = (cell_type_counts / cell_type_counts.sum()) * 100
    # Create a DataFrame
    frequencies_df = cell_type_frequencies.reset_index()
    frequencies_df.columns = ['Cell Type', 'Frequency (%)']
    frequencies_df['Group'] = group_name
    return frequencies_df

# Calculate frequencies for both groups
casp1_pos_freqs = get_frequency_percentages(mdata[mdata.obs['comparison_group'] == 'Casp1+Cx3cr1+'], 'Casp1+Cx3cr1+')
casp1_neg_freqs = get_frequency_percentages(mdata[mdata.obs['comparison_group'] == 'Casp1-Cx3cr1+'], 'Casp1-Cx3cr1+')

# Combine the two frequency dataframes
combined_freqs_df = pd.concat([casp1_pos_freqs, casp1_neg_freqs])

# Plot the grouped bar plot without enforcing an order list
plt.figure(figsize=(8, 6), dpi=1000)
sns_barplot = sns.barplot(
    x='Cell Type', 
    y='Frequency (%)', 
    hue='Group', 
    data=combined_freqs_df, 
    palette="viridis"
)

# Rotate the x labels for better readability
sns_barplot.set_xticklabels(sns_barplot.get_xticklabels(), rotation=45, ha='right')

# Set title and labels
sns_barplot.set_title('Distribution of Cx3cr1+ Immune Cells by Casp1 Expression')
sns_barplot.set_xlabel('Cell Type')
sns_barplot.set_ylabel('Frequency (%)')

# Adjust layout for better fit and disable the grid
plt.tight_layout()
plt.grid(False)
plt.savefig(os.path.join(save_dir, 'cx3cr1_casp1_distribution-v.pdf'), format='pdf', dpi=1000)  # High-resolution PNG



# DOTPLOT

sc.pl.dotplot(mdata, ['Cx3cr1', 'Casp1'], groupby = 'comparison_group', cmap = 'bwr', save = 'Cx3cr1_casp1_dotplot-bwr.pdf')

################################################################################################################################

# GENE LISTS

sc.tl.rank_genes_groups(mdata, groupby='comparison_group', method='wilcoxon', use_raw = False, layer = 'log1p')
de_results = sc.get.rank_genes_groups_df(mdata, None)
de_results = de_results[(de_results.group == 'Casp1+Cx3cr1+')].copy()

goi = pd.read_csv("/Users/maureen/Documents/projects/harris-lab/Isaac/hunter_dataset/unique-genes.csv")
goi

genes = goi['Genes'].tolist()
#genes_subset = genes[:45]
#new_elements = ['Casp1', 'Ccl5', 'Cd3g', 'Prf1', 'Cd8a']
#genes_subset = new_elements + genes_subset


de_in_goi = de_results[de_results['names'].isin(genes)]
de_in_goi['abs_logfoldchanges'] = de_in_goi['logfoldchanges'].abs()


top_50_de_goi = de_in_goi.sort_values(by=['abs_logfoldchanges', 'pvals_adj'], ascending=[False, True]).head(50)

upregulated_genes = top_50_de_goi[top_50_de_goi['logfoldchanges'] > 0]['names'].tolist()
downregulated_genes = top_50_de_goi[top_50_de_goi['logfoldchanges'] < 0]['names'].tolist()

# Concatenate gene lists to maintain group separation in the heatmap
genes_subset = upregulated_genes + downregulated_genes

# HEATMAP

#sc.pl.matrixplot(mdata, genes, groupby = 'comparison_group', dendrogram = False, use_raw = True, cmap = 'viridis', standard_scale ='var')
sc.pl.heatmap(mdata, genes_subset, groupby = 'comparison_group', dendrogram = False, use_raw = True, cmap = 'viridis', standard_scale ='var', save = 'Cx3cr1_casp1_heatmap-v.pdf')

################################################################################################################################

# MATRIX PLOT

sc.pl.matrixplot(mdata, cell_death_genes, groupby = 'comparison_group', dendrogram = True, use_raw = True, cmap = 'bwr', standard_scale = 'var')
sc.pl.matrixplot(mdata, cytokine, groupby = 'comparison_group', dendrogram = True, use_raw = True, cmap = 'bwr', standard_scale = 'var')
sc.pl.matrixplot(mdata, chemokine, groupby = 'comparison_group', dendrogram = True, use_raw = True, cmap = 'bwr', standard_scale = 'var')
sc.pl.matrixplot(mdata, cell_activation, groupby = 'comparison_group', dendrogram = True, use_raw = True, cmap = 'bwr', standard_scale = 'var') 
sc.pl.dotplot(mdata, inflammasome, groupby = 'comparison_group', dendrogram = True, use_raw = True, cmap = 'bwr')

sc.pl.dotplot(mdata, ['Cx3cr1', 'Casp1'], groupby = 'comparison_group', dendrogram = True, use_raw = True, cmap = 'bwr') 

sc.pl.matrixplot(mdata, de_genes, groupby = 'comparison_group', dendrogram = True, use_raw = True, cmap = 'bwr', standard_scale = 'var')
