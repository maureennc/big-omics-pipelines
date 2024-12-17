#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 09:23:24 2024

@author: maureen
"""


# IMPORT PACKAGES

import os
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
import seaborn as sns

import scanpy as sc
sc.set_figure_params(scanpy = True, dpi = 150, dpi_save = 400)
import squidpy as sq
import scrublet as scr

###############################################################################

# READ DATASETS

## Rivanna
#data_dir = "/scratch/mnc3ra/merfish-analysis/datasets"


## Local

data_dir = "/Users/maureen/Desktop/toxo-merfish/data/h5ad/raw"


E003 = sc.read_h5ad(os.path.join(data_dir, 'E003.h5ad'))
E007 = sc.read_h5ad(os.path.join(data_dir, 'E007.h5ad'))
E008 = sc.read_h5ad(os.path.join(data_dir, 'E008.h5ad'))
E009 = sc.read_h5ad(os.path.join(data_dir, 'E009.h5ad'))

###############################################################################

# ADDRESS BLANKS

## Check for blanks
contains_blank = E003.var_names.str.contains('Blank').any()

if contains_blank:
    print("At least one var_name contains 'Blank'")
else:
    print("No var_name contains 'Blank'")


## Filter out blanks
E003 = E003[:, ~E003.var_names.str.contains('Blank')].copy()
E007 = E007[:, ~E007.var_names.str.contains('Blank')].copy()
E008 = E008[:, ~E008.var_names.str.contains('Blank')].copy()
E009 = E009[:, ~E009.var_names.str.contains('Blank')].copy()


contains_blank = E003.var_names.str.contains('Blank').any()

if contains_blank:
    print("At least one var_name contains 'Blank'")
else:
    print("No var_name contains 'Blank'")

###############################################################################

## QC

adata_list = [E003, E007, E008, E009]

for i, adata in enumerate(adata_list):
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)


###############################################################################

# VISUALIZE QC DATA

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'volume'], palette = 'rocket')
    

###############################################################################

# FILTER CELLS BASED ON QC THRESHOLDS

## Total counts
### Floor threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts floor threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs['total_counts'] > 20, :].copy()
    print(f"Dataset {i} shape after total_counts floor threshold: {adata_list[i].shape}")


### Ceiling threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts ceiling threshold: {adata.shape}")
    total_counts_ceiling = adata.obs['total_counts'].quantile(0.97)
    adata_list[i] = adata[adata.obs['total_counts'] < total_counts_ceiling, :]
    print(f"Dataset {i} shape after total_counts ceiling threshold: {adata_list[i].shape}")

## Number of genes per cell
### Floor threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before n_genes_by_counts floor threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] > 20, :].copy()
    print(f"Dataset {i} shape after n_genes_by_counts floor threshold: {adata_list[i].shape}")

### Ceiling threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before n_genes_by_counts ceiling threshold: {adata.shape}")
    ngenes_ceiling = adata.obs['n_genes_by_counts'].quantile(0.99)
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] < ngenes_ceiling, :]
    print(f"Dataset {i} shape after n_genes_by_counts ceiling threshold: {adata_list[i].shape}")

### Volume threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before volume ceiling threshold: {adata.shape}")
    volume_ceiling = adata.obs['volume'].quantile(0.97)
    adata_list[i] = adata[adata.obs['volume'] < volume_ceiling, :]
    print(f"Dataset {i} shape after volume ceiling threshold: {adata_list[i].shape}")

for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before volume floor threshold: {adata.shape}")
    volume_floor = adata.obs['volume'].quantile(0.01)
    adata_list[i] = adata[adata.obs['volume'] > volume_floor, :]
    print(f"Dataset {i} shape after volume floor threshold: {adata_list[i].shape}")

### Anisotropy threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before anisotropy ceiling threshold: {adata.shape}")
    anisotropy_ceiling = adata.obs['anisotropy'].quantile(0.99)
    adata_list[i] = adata[adata.obs['anisotropy'] < anisotropy_ceiling, :]
    print(f"Dataset {i} shape after anisotropy ceiling threshold: {adata_list[i].shape}")

E003, E007, E008, E009 = adata_list

###############################################################################

# RECALCULATE QC METRICS
    
for i, adata in enumerate(adata_list):
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'volume'], palette = 'rocket')
    
###############################################################################

# DOUBLET DETECTION

def run_scrublet(adata_list, doublet_threshold=None):
    scrublet_rows = []
    
    for i, adata in enumerate(adata_list):
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.2)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                                  min_cells=3, 
                                                                  min_gene_variability_pctl=85, 
                                                                  n_prin_comps=30)

        if doublet_threshold is not None:
            predicted_doublets = scrub.call_doublets(threshold=doublet_threshold)
        else:
            if scrub.doublet_scores_obs_ is None or scrub.doublet_scores_sim_ is None:
                print(f"Scrublet did not return results for AnnData object {i+1}")
                continue

        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets

        # Plot Histograms
        plt.figure(figsize=(10, 6))
        sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
        sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
        plt.title(f'Scrublet Doublet Score Distribution for AnnData Object {i+1}')
        plt.xlabel('Doublet Score')
        plt.ylabel('Density')
        plt.legend()
        plt.show()

        # Extract cell barcodes
        cell_barcodes = adata.obs.index

        # Store Scrublet data with cell barcodes for each AnnData object in the list
        for barcode, obs_score, sim_score, pred_doublet in zip(cell_barcodes, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_, predicted_doublets):
            scrublet_rows.append({'Sample_Index': i+1, 
                                  'Cell_Barcode': barcode,
                                  'Observed_Score': obs_score, 
                                  'Simulated_Score': sim_score, 
                                  'Predicted_Doublet': pred_doublet})

    # Create a DataFrame from the list of rows
    scrublet_df = pd.DataFrame(scrublet_rows)
    return scrublet_df

# Run Scrublet analysis
scrublet_result = run_scrublet(adata_list, doublet_threshold=0.2)

for adata in adata_list:
    sc.pl.violin(adata, ['doublet_scores'])


###############################################################################

# CONCATENATION

## Annotate groups
E003.obs['condition'] = 'Infected'
E007.obs['condition'] = 'Infected'
E008.obs['condition'] = 'Naive'
E009.obs['condition'] = 'Naive'

## Perform concatenation
adata = E003.concatenate(E007, E008, E009, batch_key='sample', batch_categories=['E003-infected', 'E007-infected', 'E008-naive', 'E009-naive'] )
adata.obs

###############################################################################

# FINISH PRE-PROCESSING

## Normalization and Log1p transformation after concatenation
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata)
adata.layers['normalized'] = adata.X.copy()
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()

## Freeze log1p in raw slot
adata.raw = adata.copy()

## Regress & scale
sc.pp.regress_out(adata, keys='total_counts') 
sc.pp.scale(adata)
adata.layers['scaled'] = adata.X.copy()

###############################################################################


# SAVE AS H5AD

h5ad_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/filtered'
os.makedirs(h5ad_dir, exist_ok=True)

# Write the AnnData objects directly using their `write` method
E003.write(f'{h5ad_dir}/E003.h5ad')
E007.write(f'{h5ad_dir}/E007.h5ad')
E008.write(f'{h5ad_dir}/E008.h5ad')
E009.write(f'{h5ad_dir}/E009.h5ad')


adata.write(f'{h5ad_dir}/pp-concat.h5ad')
adata.obs['total_counts'].mean()
###############################################################################
