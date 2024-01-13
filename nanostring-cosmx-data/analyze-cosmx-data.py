#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze Nanostring CosMx Data

Created on Sat Jan 13 11:12:16 2024

@author: maureen
"""

# Raw data download: https://smi-public.objects.liquidweb.services/HalfBrain.zip
# Summary files download: https://smi-public.objects.liquidweb.services/Half%20%20Brain%20simple%20%20files%20.zip

##################################################################

# PACKAGES

import os
import scanpy as sc
import squidpy as sq
import scvi

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sc.set_figure_params(scanpy = True, dpi = 150, dpi_save = 400)

##################################################################

# IMPORT DATA

## Needed to create CellLabels and CompartmentLabels directories and move .tifs for each FOV into them

adata = sq.read.nanostring(path= '/project/harrislab/projects/cosmx_demo/nanostring_data/HalfBrain_20230406_205644_S1/CellStatsDir',
                           counts_file='Run1000_S1_Half_exprMat_file.csv',
                           meta_file= 'Run1000_S1_Half_metadata_file.csv',
                           fov_file= 'Run1000_S1_Half_fov_positions_file.csv'
                           )

# sc.write('demo_import.h5ad', adata)

##################################################################

# LABEL NEGATIVE CONTROL PROBES, QC

adata.var["NegPrb"] = adata.var_names.str.startswith("NegPrb")
sc.pp.calculate_qc_metrics(adata, qc_vars=["NegPrb"], inplace=True)


adata.obs["total_counts_NegPrb"].sum() / adata.obs["total_counts"].sum() * 100 # Used for FDR


## Visualize QC Metriccs

fig, axs = plt.subplots(1, 3, figsize=(15, 4))

axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)

axs[2].set_title("Transcripts per FOV")
sns.histplot(
    adata.obs.groupby("fov").sum()["total_counts"],
    kde=False,
    ax=axs[2],
)

##################################################################

# FILTERING

sc.pp.filter_genes(adata, min_cells=10)
sc.pp.filter_cells(adata, min_genes=200) # Makes a big difference
sc.pp.filter_cells(adata, min_counts=50)
sc.pp.filter_cells(adata, max_counts=7500)



##################################################################

# NORMALIZE, LOG-TRANSFORM, REGRESS, SCALE DATA

## Freeze raw data into adata.raw slot
adata.raw = adata
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, inplace=True)
adata.layers["normalized"] = adata.X.copy()

sc.pp.log1p(adata)
adata.layers["log"] = adata.X.copy()

sc.pp.regress_out(adata, ['total_counts'])
adata.layers["regressed"] = adata.X.copy()

sc.pp.scale(adata, max_value=10)
adata.layers["scaled"] = adata.X.copy()


##################################################################

# CLUSTERING (BASE SCANPY)

sc.pp.pca(adata)
sc.pp.neighbors(adata, n_pcs = 25)
sc.tl.umap(adata)
sc.tl.leiden(adata)

sc.pl.umap(adata, color = ['leiden'])

##################################################################

# CLUSTERING (scVI)

## Set up and train model
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    #categorical_covariate_keys=[''],
    continuous_covariate_keys=['total_counts'])

model = scvi.model.SCVI(adata)
model.train()

## Save model
save_dir = ('./scVI')
os.makedirs(save_dir, exist_ok=True)
model_dir = os.path.join(save_dir, "scvi_model_2024-01-11_cosmx")
model.save(model_dir, overwrite=True)

## Extract model inputs
SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

## Compute neighbors graph
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=0.3)

sc.pl.umap(
    adata,
    color=["leiden"],
    frameon=False,
)

## Perform clustering
SCVI_CLUSTERS_KEY = "leiden_scVI"
sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=0.5)

sc.pl.umap(
    adata,
    color=[SCVI_CLUSTERS_KEY],
    frameon=False,
)

## Save anndata as .h5ad
# sc.write('./h5ad/cosmx-demo-scVI-model.h5ad', adata)

##################################################################

# CELL TYPE ANNOTATION (BROAD CLASSES)

sc.tl.rank_genes_groups(adata, groupby=SCVI_CLUSTERS_KEY, method='wilcoxon', layer = 'log', use_raw = False)
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

markers_scvi = model.differential_expression(groupby = 'leiden_scVI')
markers_scvi

## Adjust and examine each cluster identifier in 'group1' (0-21)
markers_scvi[markers_scvi['group1'] == '21'].sort_values(by='lfc_mean', ascending=False).head(50).index

## Cell type dictionary
cell_type= {
"0": "Oligodendrocyte",
"1": "Neuron",
"2": "Neuron",
"3": "Neuron",
"4": "Vascular",
"5": "Neuron",
"6": "Neuron",
"7": "Neuron",
"8": "Astrocyte",
"9": "Neuron",
"10": "Astrocyte",
"11": "Neuron",
"12": "Macrophage",
"13": "OPC",
"14": "Neuron",
"15": "Neuron",
"16": "Neuron", 
"17": "Unknown",
"18": "Unknown",
"19": "Neuron",
"20": "Neuron",
"21": "Unknown"
}

adata.obs['cell_type'] = adata.obs.leiden_scVI.map(cell_type)


## Plotting
sc.pl.umap(adata, color = [SCVI_CLUSTERS_KEY], frameon=True, legend_loc='on data')

##################################################################

# CELL TYPE ANNOTATION (RESOLVING NEURONAL HETEROGENEITY)

neuron = adata[adata.obs['cell_type'] == 'Neuron'].copy()

sc.tl.rank_genes_groups(neuron, groupby=SCVI_CLUSTERS_KEY, method='wilcoxon', layer = 'log', use_raw = False)
neuron_markers = sc.get.rank_genes_groups_df(neuron, None)
neuron_markers = neuron_markers[(neuron_markers.pvals_adj < 0.05) & (neuron_markers.logfoldchanges > .5)]
neuron_markers

## Adjust and examine each neuronal cluster identifier in 'group
neuron_markers[neuron_markers['group'] == '3'].sort_values(by='logfoldchanges', ascending=False).head(50).names

## Cell class dictionary
cell_class= {
"0": "Oligodendrocyte",
"1": "Glutamatergic neuron",
"2": "GABAergic neuron",
"3": "GABAergic neuron",
"4": "Vascular",
"5": "Glutamatergic neuron",
"6": "Glutamatergic neuron",
"7": "GABAergic neuron",
"8": "Astrocyte",
"9": "GABAergic neuron",
"10": "Astrocyte",
"11": "Glutamatergic neuron",
"12": "Macrophage",
"13": "OPC",
"14": "Cholinergic neuron",
"15": "Glutamatergic neuron",
"16": "Glutamatergic neuron", 
"17": "Unknown",
"18": "Unknown",
"19": "GABAergic neuron",
"20": "Glutamatergic neuron",
"21": "Unknown"
}

adata.obs['scvi_cell_class'] = adata.obs.leiden_scVI.map(cell_class)

## Cell cluster dictionary
cell_cluster= {
"0": "Oligodendrocyte",
"1": "Glutamatergic neuron 1",
"2": "GABAergic neuron 1",
"3": "GABAergic neuron 2",
"4": "Vascular",
"5": "Glutamatergic neuron 2",
"6": "Glutamatergic neuron 3",
"7": "GABAergic neuron 3",
"8": "Astrocyte 1",
"9": "GABAergic neuron 4",
"10": "Astrocyte 2",
"11": "Glutamatergic neuron 4",
"12": "Macrophage",
"13": "OPC",
"14": "Cholinergic neuron",
"15": "Glutamatergic neuron 5",
"16": "Glutamatergic neuron 6", 
"17": "Unknown 1",
"18": "Unknown 2",
"19": "GABAergic neuron 5",
"20": "Glutamatergic neuron 7",
"21": "Unknown 3"
}

adata.obs['scvi_cell_cluster'] = adata.obs.leiden_scVI.map(cell_cluster)

## Visualization
sc.pl.umap(adata, color = 'scvi_cell_class')
sc.pl.umap(adata, color = 'scvi_cell_cluster')

##################################################################

# COMPOSITIONAL ANALYSIS

## Number of cells per cell type
cell_counts_class = adata.obs['scvi_cell_class'].value_counts().reset_index()
cell_counts_class.columns = ['scvi_cell_class', 'n_cells']
cell_counts_class = cell_counts_class.sort_values(by='n_cells', ascending=True)

plt.figure(figsize=(3, 3))
sns.barplot(data=cell_counts_class, x='scvi_cell_class', y='n_cells', palette='rocket')
plt.xticks(rotation=90)
plt.xlabel('Cell type')
plt.ylabel('Number of Cells')
plt.title('Number of Cells per Cell Type')
plt.grid(False)
plt.show()

## Frequency of cells per cell type

## Bar plot
total_cells = cell_counts_class['n_cells'].sum()
cell_counts_class['frequency'] = (cell_counts_class['n_cells'] / total_cells) * 100

plt.figure(figsize=(3, 3))
sns.barplot(data=cell_counts_class, x='scvi_cell_class', y='frequency', palette='rocket')
plt.xticks(rotation=90)
plt.xlabel('Cell type')
plt.ylabel('Frequency (%)')
plt.title('Frequency of Cells per Cell Type')
plt.grid(False)
plt.show()

### Pie Chart
plt.figure(figsize=(8, 8))
pie = plt.pie(cell_counts_class['n_cells'], autopct='%1.1f%%', startangle=140)

# Create the legend
plt.legend(pie[0], cell_counts_class['scvi_cell_class'], title="Cell type", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
plt.title('Distribution of Cells per Cell Type')
plt.show()


## Number of cells per cluster
cell_counts_cluster = adata.obs['scvi_cell_cluster'].value_counts().reset_index()
cell_counts_cluster.columns = ['scvi_cell_cluster', 'n_cells']

plt.figure(figsize=(5, 3))
sns.barplot(data=cell_counts_cluster, x='scvi_cell_cluster', y='n_cells', palette='rocket')
plt.xticks(rotation=90)
plt.xlabel('Cluster')
plt.ylabel('Number of Cells')
plt.title('Number of Cells per Cluster')
plt.grid(False)
plt.show()

##################################################################

# DIFFERENTIAL EXPRESSION (DOTPLOT)

## Toggle through groupings to find good marker genes for different subsets of neurons
sc.tl.rank_genes_groups(adata, groupby='scvi_cell_class', groups=['Cholinergic neuron'], reference='Glutamatergic neuron', use_raw = False, layer = 'log')
sc.pl.rank_genes_groups(adata, n_genes=20)


## Plotting
dotplot_genes = ['Plp1', 'Mag', 'Mog', 'Olig1', # Oligodendrocytes
                 'Gja1', 'Apoe', 'Atp1a2', 'Aqp4', # Astrocytes
                 'Pdgfra', # OPCs
                 'Hexb', 'Csf1r', 'C1qa', # Macrophages
                 'Cldn5', 'Vtn', 'Rgs5', 'Flt1', # Vascular
                 'Rbfox3', 'Snap25', 'Thy1', # Pan-neuronal
                 'Chat', # Cholinergic neurons
                 'Gad1', 'Gad2', # GABAergic neurons
                 'Slc17a7', 'Mef2c', 'Camk2a'] # Glutamatergic neurons

sc.pl.dotplot(adata, dotplot_genes, groupby = 'scvi_cell_class', dendrogram = True, layer = 'scaled')

##################################################################

# EXPORT

# sc.write('./h5ad/cosmx-demo-clustered.h5ad', adata)

##################################################################