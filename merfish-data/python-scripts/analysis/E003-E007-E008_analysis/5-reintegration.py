#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 21:16:06 2024

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

data_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/h5ad/processed-data'

adata_concat = sc.read_h5ad(os.path.join(data_dir, 'concat-model_A-roi-annotated.h5ad'))
adata_infected = sc.read_h5ad(os.path.join(data_dir, 'infected-model_A-roi-annotated.h5ad'))
#adata_naive = sc.read_h5ad(os.path.join(data_dir, 'naive-model_A-roi-annotated.h5ad'))

################################################################################################################################

# SET UP AND TRAIN MODEL B, CONCAT

## Model B
scvi.model.SCVI.setup_anndata(
    adata_concat,
    layer='counts',
    categorical_covariate_keys=['sample', 'condition'],
    continuous_covariate_keys=['total_counts', 'doublet_scores', 'centroid_theta', 'centroid_r', 'anisotropy'],
)

model_B = scvi.model.SCVI(adata_concat)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_B.train()

################################################################################################################################

# SET UP AND TRAIN MODEL C, INFECTED

## Model C
scvi.model.SCVI.setup_anndata(
    adata_infected,
    layer='counts',
    categorical_covariate_keys=['sample', 'condition'],
    continuous_covariate_keys=['total_counts', 'doublet_scores', 'centroid_theta', 'centroid_r', 'anisotropy'],
)

model_C = scvi.model.SCVI(adata_infected)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_C.train()


################################################################################################################################

# SAVE / IMPORT MODEL

scvi_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/scvi'

## Save model_B
#model_B_dir = os.path.join(scvi_dir, 'model_B')
#print(model_B_dir)
#model_B.save(model_B_dir)

## Save model_C
#model_C_dir = os.path.join(scvi_dir, 'model_C')
#print(model_C_dir)
#model_C.save(model_C_dir)


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

# EVALUATE TRAINED MODELS

## model B
training_history = model_B.history
training_history

training_history_df = pd.DataFrame(index=training_history['kl_weight'].index)

for key, df in training_history.items():
    training_history_df = training_history_df.join(df, how='outer')

## Visualize results
training_history_df.reset_index(inplace=True)

plt.figure(figsize=(5, 20))
## ELBO
plt.subplot(3, 1, 1)
plt.plot(training_history_df['epoch'], training_history_df['elbo_train'], label='ELBO')
plt.xlabel('Epochs')
plt.ylabel('ELBO')
plt.title('ELBO over Training Epochs')
plt.legend()

## Training Loss
plt.subplot(3, 1, 2)
plt.plot(training_history_df['epoch'], training_history_df['train_loss_epoch'], label='Training Loss')
plt.xlabel('Epochs')
plt.ylabel('Training Loss')
plt.title('Training Loss over Epochs')
plt.legend()

## KL Divergence (Local)
plt.subplot(3, 1, 3)
plt.plot(training_history_df['epoch'], training_history_df['kl_local_train'], label='KL Divergence (Local)')
plt.xlabel('Epochs')
plt.ylabel('KL Divergence (Local)')
plt.title('KL Divergence over Epochs')
plt.legend()

## Adjust layout
plt.tight_layout()
plt.show()



## Model C
training_history = model_C.history
training_history

training_history_df = pd.DataFrame(index=training_history['kl_weight'].index)

for key, df in training_history.items():
    training_history_df = training_history_df.join(df, how='outer')

## Visualize results
training_history_df.reset_index(inplace=True)

plt.figure(figsize=(5, 20))
## ELBO
plt.subplot(3, 1, 1)
plt.plot(training_history_df['epoch'], training_history_df['elbo_train'], label='ELBO')
plt.xlabel('Epochs')
plt.ylabel('ELBO')
plt.title('ELBO over Training Epochs')
plt.legend()

## Training Loss
plt.subplot(3, 1, 2)
plt.plot(training_history_df['epoch'], training_history_df['train_loss_epoch'], label='Training Loss')
plt.xlabel('Epochs')
plt.ylabel('Training Loss')
plt.title('Training Loss over Epochs')
plt.legend()

## KL Divergence (Local)
plt.subplot(3, 1, 3)
plt.plot(training_history_df['epoch'], training_history_df['kl_local_train'], label='KL Divergence (Local)')
plt.xlabel('Epochs')
plt.ylabel('KL Divergence (Local)')
plt.title('KL Divergence over Epochs')
plt.legend()

## Adjust layout
plt.tight_layout()
plt.show()


################################################################################################################################

# EXTRACT LATENT REPRESENTATION

## Model B
SCVI_LATENT_KEY = "X_scVI"
latent = model_B.get_latent_representation()
adata_concat.obsm[SCVI_LATENT_KEY] = latent
latent.shape

adata_concat.layers['scvi_normalized'] = model_B.get_normalized_expression(library_size = 1e4)
adata_concat.layers

## Model C
SCVI_LATENT_KEY = "X_scVI"
latent = model_C.get_latent_representation()
adata_infected.obsm[SCVI_LATENT_KEY] = latent
latent.shape

adata_infected.layers['scvi_normalized'] = model_C.get_normalized_expression(library_size = 1e4)
adata_infected.layers

################################################################################################################################

# INITIAL CLUSTERING

## Concat
sc.pp.neighbors(adata_concat, use_rep = 'X_scVI', random_state = 0) # Use latent representation to build neighbors graph
sc.tl.umap(adata_concat, min_dist = 0.3)
sc.tl.leiden(adata_concat, key_added='leiden_scVI', resolution=1)

sc.pl.umap(adata_concat, color = ['leiden_scVI'], legend_loc = 'on data')
sc.pl.umap(adata_concat, color = ['leiden_scVI', 'condition', 'sample'])
sc.pl.umap(adata_concat, color = ['P2ry12', 'Tnfaip2', 'Ccr2', 'Cd4', 'Cd8a', 'Rbfox3', 'Olig1', 'Pdgfra', 'Rgs5', 'Pecam1', 'Tie1', 'Aqp4', 'Slc17a6', 'Slc17a7', 'Gad2', 'Gas6'])


## Infected
sc.pp.neighbors(adata_infected, use_rep = 'X_scVI', random_state = 0) # Use latent representation to build neighbors graph
sc.tl.umap(adata_infected, min_dist = 0.3)
sc.tl.leiden(adata_infected, key_added='leiden_scVI', resolution=1)

sc.pl.umap(adata_infected, color = ['leiden_scVI'], legend_loc = 'on data')
sc.pl.umap(adata_infected, color = ['leiden_scVI', 'condition', 'sample'])
sc.pl.umap(adata_infected, color = ['P2ry12', 'Tnfaip2', 'Ccr2', 'Cd4', 'Cd8a', 'Rbfox3', 'Olig1', 'Pdgfra', 'Rgs5', 'Pecam1', 'Tie1', 'Aqp4', 'Slc17a6', 'Slc17a7', 'Gad2', 'Gas6'])

################################################################################################################################

# EXPORT

save_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/h5ad/processed-data'

adata_concat.write_h5ad(os.path.join(save_dir, 'concat-model_B.h5ad'))
adata_infected.write_h5ad(os.path.join(save_dir, 'infected-only-model_B.h5ad'))

################################################################################################################################
