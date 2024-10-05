#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 12:14:54 2024

@author: maureen
"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
import random
import torch
import scvi

################################################################################

# SETTINGS

## Random seed
random.seed(0)
torch.manual_seed(0)
np.random.seed(0)
scvi.settings.seed = 0

###############################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/3/h5ad/processed/"

adata = sc.read_h5ad(os.path.join(data_dir, '1-syk-pp-concat-v3.h5ad'))

###############################################################################

# ASSESS COVARIATE MULTILINEARITY

obs_df = adata.obs.copy()

correlation_matrix = obs_df[['total_counts', 'n_genes_by_counts', 'cell_ID', 'fov', 'doublet_scores', 'pct_counts_NegPrb']].corr()
cm = pd.DataFrame(correlation_matrix)

###############################################################################

# SCVI TRAINING

scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    categorical_covariate_keys=['orig_dataset'],
    continuous_covariate_keys=['total_counts'],
)

model_A = scvi.model.SCVI(adata)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_A.train()

###############################################################################

# SAVE / IMPORT MODEL

scvi_dir = "/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/3/scvi"

## Save 
model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
#model_A.save(model_A_dir)


## Import model
scvi_dir = "/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/3/scvi"

model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
model_A = scvi.model.SCVI.load(model_A_dir, adata=adata)
model_A

###############################################################################

# EVALUATE TRAINED MODEL A

training_history = model_A.history
training_history


training_history_df = pd.DataFrame(index=training_history['kl_weight'].index)

for key, df in training_history.items():
    training_history_df = training_history_df.join(df, how='outer')

training_history_df.reset_index(inplace=True)

plt.figure(figsize=(5, 20))

plt.subplot(3, 1, 1)
plt.plot(training_history_df['epoch'], training_history_df['elbo_train'], label='ELBO')
plt.xlabel('Epochs')
plt.ylabel('ELBO')
plt.title('ELBO over Training Epochs')
plt.legend()


plt.subplot(3, 1, 2)
plt.plot(training_history_df['epoch'], training_history_df['train_loss_epoch'], label='Training Loss')
plt.xlabel('Epochs')
plt.ylabel('Training Loss')
plt.title('Training Loss over Epochs')
plt.legend()


plt.subplot(3, 1, 3)
plt.plot(training_history_df['epoch'], training_history_df['kl_local_train'], label='KL Divergence (Local)')
plt.xlabel('Epochs')
plt.ylabel('KL Divergence (Local)')
plt.title('KL Divergence over Epochs')
plt.legend()
plt.tight_layout()
plt.show()

###############################################################################

# EXTRACT LATENT REPRESENTATION

adata.obsm
SCVI_LATENT_KEY = "X_scVI"

latent_A = model_A.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent_A
latent_A.shape

adata.layers['scvi_normalized'] = model_A.get_normalized_expression()
adata.layers

###############################################################################

# INITIAL CLUSTERING

sc.pp.neighbors(adata, use_rep = 'X_scVI', random_state = 0) 
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added='leiden_scVI', resolution=1)

sc.pl.umap(adata, color = ['leiden_scVI'], legend_loc = 'on data')


sc.tl.rank_genes_groups(adata, groupby='leiden_scVI', method='wilcoxon', use_raw = False, layer = 'log1p')
sc.pl.rank_genes_groups(adata, n_genes = 20, fontsize = 20)


genes = ['Nrgn', 'Rbfox3', 'Snap25', 'Olfm1', 'Camk2a', 'Camk2b', 'Calm1', 'Calm3', 'Gria2', 
         'Tmsb4x', 'Gad1', 'Thy1', 'Pcp4',
         'Cd44', 'Tgfbr1', 'Itgax', 'Cd3d', 'Rgs5', 
         'Apoe', 'Slc1a2', 'Slc1a3', 'Gfap', 'Glul',
         'Mbp', 'Plp1',
         'Nap1l5', 'Aldoa', 'Bsg',
         'Cst3', 'Hexb', 'P2ry12', 'Csf1r', 'Olig1', 'Olig2', 'Pdgfra', 'Cspg5', 
        'Nnat', 'Rgs5', 'Pecam1']

sc.tl.dendrogram(adata, groupby = 'leiden_scVI')
sc.pl.dotplot(adata, genes, groupby = 'leiden_scVI', standard_scale = 'var', use_raw = True, dendrogram = True)

sc.pl.umap(adata, color = ['leiden_scVI', 'Nrgn', 'Rbfox3', 'Snap25', 'Gad1', 'Mbp', 'Plp1', 'Hexb', 'P2ry12', 'Pdgfra', 'Olig1', 'Aqp4'], legend_loc = 'on data')

###############################################################################

# CELL TYPING

cell_type= { 
"0": "Neuron",
"1": "Astrocyte",
"2": "Neuron",
"3": "Neuron",
"4": "Oligodendrocyte",
"5": "Unknown",
"6": "Neuron",
"7": "Vascular",
"8": "Unknown",
"9": "Unknown",
"10": "Microglia",
"11": "OPC"
}

adata.obs['cell_type'] = adata.obs.leiden_scVI.map(cell_type)


cluster= { 
"0": "Neuron 1",
"1": "Astrocyte",
"2": "Neuron 2",
"3": "Neuron 3",
"4": "Oligodendrocyte",
"5": "Unknown, neuronal 1",
"6": "Neuron 4",
"7": "Vascular",
"8": "Unknown, neuronal 2",
"9": "Unknown, glial",
"10": "Microglia",
"11": "OPC"
}

adata.obs['cluster'] = adata.obs.leiden_scVI.map(cluster)

###############################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/3/h5ad/"

adata.X = csr_matrix(adata.X)
adata.write_h5ad(os.path.join(save_dir, '2-syk-integrated-annotated.h5ad'))

###############################################################################