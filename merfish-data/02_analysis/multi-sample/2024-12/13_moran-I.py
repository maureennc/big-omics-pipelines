#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 18:39:36 2024

@author: maureen
"""

import squidpy as sq
import scanpy as sc
import os
import matplotlib.pyplot as plt
import scipy


###############################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/processed'

adata = sc.read_h5ad(os.path.join(data_dir, 'concat-split-trained-annotated.h5ad'))

###############################################################################

# MORAN I

E003 = adata[adata.obs['sample'] == 'E003-infected'].copy()
E007 = adata[adata.obs['sample'] == 'E007-infected'].copy()
E008 = adata[adata.obs['sample'] == 'E008-naive'].copy()
E009 = adata[adata.obs['sample'] == 'E009-naive'].copy()

adata_list = [E003, E007, E008, E009]


sample_condition_map = {
    'E003': 'Infected',
    'E007': 'Infected',
    'E008': 'Naive',
    'E009': 'Naive'
}

adata_list = [E003, E007, E008, E009]

all_moran_results = []

for i, anndata in enumerate(adata_list):
    sq.gr.spatial_neighbors(anndata)
    sq.gr.spatial_autocorr(anndata, mode='moran', n_perms=100, n_jobs=1)
    moran_df = anndata.uns['moranI'].copy()
    moran_df.reset_index(inplace=True)  # Reset index to preserve gene names
    moran_df['sample_id'] = f'Sample_{i + 1}'
    all_moran_results.append(moran_df)

combined_moran_df = pd.concat(all_moran_results, ignore_index=True)



## Update DF with annotations

sample_condition_map = {
    0: ('E003', 'Infected'),
    1: ('E007', 'Infected'),
    2: ('E008', 'Naive'),
    3: ('E009', 'Naive')
}

combined_moran_df['sample_index'] = combined_moran_df['sample_id'].str.extract(r'(\d+)').astype(int) - 1
combined_moran_df[['sample_name', 'condition']] = combined_moran_df['sample_index'].map(sample_condition_map).apply(pd.Series)
combined_moran_df.drop(columns=['sample_index'], inplace=True)

print(combined_moran_df.head())

#combined_moran_df.to_csv('moranI_results.csv', index=False)


###############################################################################

# SCATTERPLOT - MORAN I VS P VALUE


# Calculate mean Moran's I and p-values per gene
mean_values = combined_moran_df.groupby('index').agg(
    mean_I=('I', 'mean'),
    mean_pval=('pval_norm', 'mean')
).reset_index()

# Scatterplot of Moran's I vs. p-values
plt.figure(figsize=(8, 6))
sns.scatterplot(data=mean_values, x='mean_I', y='mean_pval', hue='index', legend=False, s=100)
plt.yscale('log')  # Log-scale for p-values
plt.xlabel("Mean Moran's I")
plt.ylabel("Mean p-value (log scale)")
plt.title("Moran's I vs. p-value")
plt.show()

###############################################################################

# HEATMAP - TOP MORAN I ACROSS SAMPLES

plt.rcParams['figure.dpi'] = 500
plt.rcParams['font.size'] = 12

top_genes = combined_moran_df.groupby('index')['I'].mean().nlargest(50).index
top_genes_df = combined_moran_df[combined_moran_df['index'].isin(top_genes)]

# Pivot data to create a heatmap with individual samples
heatmap_data = top_genes_df.pivot_table(values='I', index='index', columns='sample_name')

# Sort the heatmap data by mean Moran's I in descending order
heatmap_data = heatmap_data.loc[heatmap_data.mean(axis=1).sort_values(ascending=False).index]

ordered_columns = ['E008', 'E009', 'E003', 'E007']
heatmap_data = heatmap_data[ordered_columns]

# Plot heatmap
plt.figure(figsize=(10, 12))
sns.heatmap(
    heatmap_data, annot=True, cmap='rocket', cbar_kws={'label': 'Moran’s I'}
)
plt.title('')
plt.ylabel(' ')
plt.xlabel('')

# Rotate the y-axis labels to horizontal
plt.yticks(rotation=0)
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.show()

###############################################################################

# BAR PLOT - DELTA

# Filter Naive and Infected samples separately
naive_df = combined_moran_df[combined_moran_df['condition'] == 'Naive']
infected_df = combined_moran_df[combined_moran_df['condition'] == 'Infected']

# Calculate the mean Moran's I for each gene across all samples
naive_means = naive_df.groupby('index')['I'].mean().reset_index(name='naive_I')
infected_means = infected_df.groupby('index')['I'].mean().reset_index(name='infected_I')

# Merge the Naive and Infected data on the 'index' (gene name)
delta_df = pd.merge(naive_means, infected_means, on='index', how='inner')

# Calculate the delta (Infected - Naive)
delta_df['delta_I'] = delta_df['infected_I'] - delta_df['naive_I']

# Debug: Check if Cxcl10 is present
print(delta_df[delta_df['index'] == 'Cxcl10'])

# Sort the delta values in descending order and select the top 10
top_delta_genes = delta_df.sort_values(by='delta_I', ascending=False).head(10)

# Plot the top 10 genes with the largest delta
plt.figure(figsize=(10, 6))
top_delta_genes.set_index('index')['delta_I'].plot(kind='bar', color='coral')
plt.title('Top 10 Genes with Largest Moran’s I Changes (Infected - Naive)')
plt.ylabel('Delta Moran’s I')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()


###############################################################################

# BAR CHART - GROUPED

grouped_data = combined_moran_df.pivot_table(
    values='I', index='index', columns='sample_name'
)

# Calculate delta values and add to the grouped data
grouped_data['delta_I'] = (
    grouped_data[['E003', 'E007']].mean(axis=1) -
    grouped_data[['E008', 'E009']].mean(axis=1)
)

# Select top 10 genes with the largest delta values
top_genes = grouped_data['delta_I'].nlargest(50).index
top_data = grouped_data.loc[top_genes].drop(columns='delta_I')

# Extract two colors from the 'rocket' colormap
rocket_colors = sns.color_palette('rocket', 2)
colors = [rocket_colors[1], rocket_colors[1], rocket_colors[0], rocket_colors[0]]  # Naive first

# Plot grouped bar chart
fig, ax = plt.subplots(figsize=(10, 8))
top_data.plot(kind='bar', ax=ax, color=colors)

plt.title('')
plt.ylabel('Moran’s I')
plt.xticks(rotation=45, ha='right')

# Customize legend inside the plot
handles = [
    plt.Line2D([0], [0], marker='o', color=rocket_colors[1], markersize=10, linestyle='', label='Naive'),
    plt.Line2D([0], [0], marker='o', color=rocket_colors[0], markersize=10, linestyle='', label='Infected')
]
plt.legend(
    handles=handles, title='Condition', 
    loc='upper right', bbox_to_anchor=(0.98, 0.95), frameon=True
)

plt.tight_layout()
plt.show()



###############################################################################


# SPATIAL SCATTERPLOTS

plt.rcParams['figure.dpi'] = 200


sq.gr.spatial_neighbors(E003)
sq.gr.spatial_neighbors(E007)
sq.gr.spatial_neighbors(E008)
sq.gr.spatial_neighbors(E009)


# Cxcl10
sq.pl.spatial_scatter(E003,color = 'Cxcl10', shape=None, cmap = 'rocket', vmax = 4)
sq.pl.spatial_scatter(E007,color = 'Cxcl10', shape=None, cmap = 'rocket', vmax = 4)
sq.pl.spatial_scatter(E008,color = 'Cxcl10', shape=None, cmap = 'rocket', vmax = 4)
sq.pl.spatial_scatter(E009,color = 'Cxcl10', shape=None, cmap = 'rocket', vmax = 4)

# Mef2c
sq.pl.spatial_scatter(E003,color = 'Mef2c', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E007,color = 'Mef2c', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E008,color = 'Mef2c', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E009,color = 'Mef2c', shape=None, cmap = 'rocket', frameon = False)

# Nr4a2
sq.pl.spatial_scatter(E003,color = 'Nr4a2', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E007,color = 'Nr4a2', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E008,color = 'Nr4a2', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E009,color = 'Nr4a2', shape=None, cmap = 'rocket', frameon = False)

# Il33
sq.pl.spatial_scatter(E003,color = 'Il33', shape=None, cmap = 'rocket', frameon = False, vmax = 4)
sq.pl.spatial_scatter(E007,color = 'Il33', shape=None, cmap = 'rocket', frameon = False, vmax = 4)
sq.pl.spatial_scatter(E008,color = 'Il33', shape=None, cmap = 'rocket', frameon = False, vmax = 4)
sq.pl.spatial_scatter(E009,color = 'Il33', shape=None, cmap = 'rocket', frameon = False, vmax = 4)

# Stat1
sq.pl.spatial_scatter(E003,color = 'Cybb', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E007,color = 'Cybb', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E008,color = 'Cybb', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E009,color = 'Stat1', shape=None, cmap = 'rocket', frameon = False)

# Stat1
sq.pl.spatial_scatter(E003,color = 'Ptprc', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E007,color = 'Ptprc', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E008,color = 'Ptprc', shape=None, cmap = 'rocket', frameon = False)
sq.pl.spatial_scatter(E009,color = 'Ptprc', shape=None, cmap = 'rocket', frameon = False)


# Cxcl10
sq.pl.spatial_scatter(E003, color='Cxcl10', shape=None, cmap='rocket', frameon=False, vmax=4)
sq.pl.spatial_scatter(E007, color='Cxcl10', shape=None, cmap='rocket', frameon=False, vmax=4)
sq.pl.spatial_scatter(E008, color='Cxcl10', shape=None, cmap='rocket', frameon=False, vmax=4)
sq.pl.spatial_scatter(E009, color='Cxcl10', shape=None, cmap='rocket', frameon=False, vmax=4)

# Stat1
sq.pl.spatial_scatter(E003, color='Stat1', shape=None, cmap='rocket', frameon=False, vmax=4)
sq.pl.spatial_scatter(E007, color='Stat1', shape=None, cmap='rocket', frameon=False, vmax=4)
sq.pl.spatial_scatter(E008, color='Stat1', shape=None, cmap='rocket', frameon=False, vmax=4)
sq.pl.spatial_scatter(E009, color='Stat1', shape=None, cmap='rocket', frameon=False, vmax=4)

# Ccl2
sq.pl.spatial_scatter(E003, color='Ccl2', shape=None, cmap='rocket', frameon=False, vmax=3)
sq.pl.spatial_scatter(E007, color='Ccl2', shape=None, cmap='rocket', frameon=False, vmax=3)
sq.pl.spatial_scatter(E008, color='Ccl2', shape=None, cmap='rocket', frameon=False, vmax=3)
sq.pl.spatial_scatter(E009, color='Ccl2', shape=None, cmap='rocket', frameon=False, vmax=3)

# C3
sq.pl.spatial_scatter(E003, color='C3', shape=None, cmap='rocket', frameon=False, vmax=4)
sq.pl.spatial_scatter(E007, color='C3', shape=None, cmap='rocket', frameon=False, vmax=4)
sq.pl.spatial_scatter(E008, color='C3', shape=None, cmap='rocket', frameon=False, vmax=4)
sq.pl.spatial_scatter(E009, color='C3', shape=None, cmap='rocket', frameon=False, vmax=4)

# Ptprc
sq.pl.spatial_scatter(E003, color='Ptprc', shape=None, cmap='rocket', frameon=False, vmax=3)
sq.pl.spatial_scatter(E007, color='Ptprc', shape=None, cmap='rocket', frameon=False, vmax=3)
sq.pl.spatial_scatter(E008, color='Ptprc', shape=None, cmap='rocket', frameon=False, vmax=3)
sq.pl.spatial_scatter(E009, color='Ptprc', shape=None, cmap='rocket', frameon=False, vmax=3)

# List of genes and samples to plot
genes = ['Cxcl10', 'Stat1', 'Nos2', 'C3', 'Ptprc', 'Ccl2']
samples = [E003, E007, E008, E009]

# Loop through each gene and sample to generate scatter plots
for gene in genes:
    for sample in samples:
        sq.pl.spatial_scatter(
            sample, color=gene, shape=None, cmap='rocket', frameon=False
        , vmax = 3)



###############################################################################

