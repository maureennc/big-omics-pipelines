#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:58:11 2024

@author: maureen
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import scrublet as scr
from scipy.sparse import csr_matrix

###############################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/3/h5ad/raw/"

syk2 = sc.read_h5ad(os.path.join(data_dir, "syk2-v2-raw.h5ad"))
syk3 = sc.read_h5ad(os.path.join(data_dir, "syk3-v2-raw.h5ad"))

annotations_path = '/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/3/annotation/csv'
syk2_anno = pd.read_csv(os.path.join(annotations_path, 'syk2-annotations.csv'))
syk3_anno = pd.read_csv(os.path.join(annotations_path, 'syk3-annotations.csv'))

###############################################################################

# ANNOTATIONS

syk2_anno.set_index('fov_anno', inplace=True)
syk2.obs['fov'] = syk2.obs['fov'].astype(str)
syk2.obs['fov_anno'] = 'FOV' + syk2.obs['fov'].astype(str).str.zfill(3)
syk2.obs = syk2.obs.join(syk2_anno[['dataset_anno', 'sample_number', 'region', 'group']], on='fov_anno')



syk3_anno.set_index('fov_anno', inplace=True)
syk3.obs['fov'] = syk3.obs['fov'].astype(str)
syk3.obs['fov_anno'] = 'FOV' + syk3.obs['fov'].astype(str).str.zfill(3)
syk3.obs = syk3.obs.join(syk3_anno[['dataset_anno', 'sample_number', 'region', 'group']], on='fov_anno')

condition= { 
"A": "Mock, WT",
"B": "Mock, Syk-KO",
"C": "SARS-CoV-2, WT",
"D": "SARS-CoV-2, Syk-KO",
}

syk2.obs['condition'] = syk2.obs.group.map(condition)

order = ["Mock, WT",
         "Mock, Syk-KO",
         "SARS-CoV-2, WT",
         "SARS-CoV-2, Syk-KO",
         "Excluded"
]
syk2.obs['condition'] = pd.Categorical(syk2.obs['condition'], categories=order, ordered=True)


syk3.obs['condition'] = syk3.obs.group.map(condition)

order = ["Mock, WT",
         "Mock, Syk-KO",
         "SARS-CoV-2, WT",
         "SARS-CoV-2, Syk-KO",
         "Excluded"
]
syk3.obs['condition'] = pd.Categorical(syk3.obs['condition'], categories=order, ordered=True)

###############################################################################

# QC

adata_list = [syk2, syk3]

for i, adata in enumerate(adata_list):
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    
# LABEL NEGATIVE CONTROL PROBES, QC

syk2.var["NegPrb"] = syk2.var_names.str.startswith("Negative")
sc.pp.calculate_qc_metrics(syk2, qc_vars=["NegPrb"], inplace=True)
syk2.obs["total_counts_NegPrb"].sum() / syk2.obs["total_counts"].sum() * 100 # Used for FDR


syk3.var["NegPrb"] = syk3.var_names.str.startswith("Negative")
sc.pp.calculate_qc_metrics(syk3, qc_vars=["NegPrb"], inplace=True)
syk3.obs["total_counts_NegPrb"].sum() / syk3.obs["total_counts"].sum() * 100 # Used for FDR

negative_mask = syk2.var.index.str.startswith('Negative')
system_control_mask = syk2.var.index.str.startswith('SystemControl')
combined_mask = negative_mask | system_control_mask
syk2 = syk2[:, ~combined_mask].copy()

negative_mask = syk3.var.index.str.startswith('Negative')
system_control_mask = syk3.var.index.str.startswith('SystemControl')
combined_mask = negative_mask | system_control_mask
syk3 = syk3[:, ~combined_mask].copy()

adata_list = [syk2, syk3]

###############################################################################

# VISUALIZE QC DATA

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts'])
    

###############################################################################

# FILTER CELLS BASED ON QC THRESHOLDS

for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts floor threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs['total_counts'] > 100, :].copy()
    print(f"Dataset {i} shape after total_counts floor threshold: {adata_list[i].shape}")

for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts ceiling threshold: {adata.shape}")
    total_counts_ceiling = adata.obs['total_counts'].quantile(0.95)
    adata_list[i] = adata[adata.obs['total_counts'] < total_counts_ceiling, :]
    print(f"Dataset {i} shape after total_counts ceiling threshold: {adata_list[i].shape}")

for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before n_genes_by_counts ceiling threshold: {adata.shape}")
    ngenes_ceiling = adata.obs['n_genes_by_counts'].quantile(0.975)
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] < ngenes_ceiling, :]
    print(f"Dataset {i} shape after n_genes_by_counts ceiling threshold: {adata_list[i].shape}")


syk2, syk3 = adata_list

###############################################################################

# RECALCULATE QC METRICS

for i, adata in enumerate(adata_list):
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

qc_metrics = ['total_counts', 'n_genes_by_counts']

for i, adata in enumerate(adata_list, start=1):
    print(f"Dataset {i}:")
    for metric in qc_metrics:
        mean_value = adata.obs[metric].mean()
        print(f"Mean {metric}: {mean_value}")
    print("-" * 30) 

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts'])


###############################################################################

# DOUBLET DETECTION

def run_scrublet(adata_list):

    scrublet_rows = []

    for i, adata in enumerate(adata_list):

        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.10, random_state = 0)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                                  min_cells=3, 
                                                                  min_gene_variability_pctl=85, 
                                                                  n_prin_comps=50)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets

        plt.figure(figsize=(10, 6))
        sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
        sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
        plt.title(f'Scrublet Doublet Score Distribution for Sample {i+1}')
        plt.xlabel('Doublet Score')
        plt.ylabel('Density')
        plt.legend()
        plt.grid(False)
        plt.show()

        cell_barcodes = adata.obs.index

        for barcode, obs_score, sim_score, pred_doublet in zip(cell_barcodes, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_, predicted_doublets):
            scrublet_rows.append({'Sample_Index': i+1, 
                                  'Cell_Barcode': barcode,
                                  'Observed_Score': obs_score, 
                                  'Simulated_Score': sim_score, 
                                  'Predicted_Doublet': pred_doublet})

    scrublet_df = pd.DataFrame(scrublet_rows)
    return scrublet_df

scrublet_result = run_scrublet(adata_list)


###############################################################################

# FILTER DOUBLETS

def filter_doublets(adata_list):
    filtered_list = []
    for i, adata in enumerate(adata_list):
        print(f"Dataset {i + 1} shape before doublet filter: {adata.shape}")
        filtered_adata = adata[adata.obs['predicted_doublets'] == False].copy()
        filtered_list.append(filtered_adata)
    return filtered_list

adata_list = filter_doublets(adata_list)

syk2, syk3 = adata_list

for i, adata in enumerate(adata_list):
    print(f"Dataset {i + 1} shape after doublet filter: {adata.shape}")
    

###############################################################################

# FINISH PRE-PROCESSING

adata = syk2.concatenate(syk3, batch_key='orig_dataset', batch_categories=['syk2', 'syk3'])

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata)
adata.layers['normalized'] = adata.X.copy()
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()
adata.raw = adata

## Scale data
sc.pp.scale(adata, max_value=10)
adata.layers['scaled'] = adata.X.copy()

#sc.pp.neighbors(adata, random_state = 0)
#sc.tl.umap(adata)
#sc.tl.leiden(adata, resolution=1)
#sc.pl.umap(adata, color = ['leiden'])

###############################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/3/h5ad/processed"

adata.X = csr_matrix(adata.X)
adata.write_h5ad(os.path.join(save_dir, '1-syk-pp-concat-v3.h5ad'))

###############################################################################