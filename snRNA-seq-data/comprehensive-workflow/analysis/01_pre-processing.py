#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:40:38 2024

@author: maureen
"""

import os
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'

import seaborn as sns
import pandas as pd
import scanpy as sc
#sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

from scipy.sparse import csr_matrix
import scrublet as scr

# environment = 'sc-pp'

###############################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/azenta-files/01_analysis/cellranger_count"

sample_A = sc.read_10x_mtx(os.path.join(data_dir, "1/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_B = sc.read_10x_mtx(os.path.join(data_dir, "3/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_C = sc.read_10x_mtx(os.path.join(data_dir, "5/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_D = sc.read_10x_mtx(os.path.join(data_dir, "6/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_E = sc.read_10x_mtx(os.path.join(data_dir, "7/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_F = sc.read_10x_mtx(os.path.join(data_dir, "8/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)

###############################################################################

## CALCULATE MITOCHONDRIAL AND RIBOSOMAL GENES

adata_list = [sample_A, sample_B, sample_C, sample_D, sample_E, sample_F]

for i, adata in enumerate(adata_list):
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribosomal'] = adata.var_names.str.match('^(Rpl|Rps)\\d+')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)

###############################################################################

# VISUALIZE QC DATA

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal'])
    

os.chdir('/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/qc/raw')

## Map each AnnData object to the specific group identifier
identifiers = ['A', 'B', 'C', 'D', 'E', 'F']
adata_list = [sample_A, sample_B, sample_C, sample_D, sample_E, sample_F]

## Save plots
for adata, id in zip(adata_list, identifiers):
    keys = ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal']
    for key in keys:
        sc.pl.violin(adata, key, show=False)
        plt.savefig(f"{key}_violin_group_{id}_raw.pdf")
        plt.close()

###############################################################################

# FILTER CELLS BASED ON QC THRESHOLDS

## Total counts
### Floor threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts floor threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs['total_counts'] > 1000, :].copy()
    print(f"Dataset {i} shape after total_counts floor threshold: {adata_list[i].shape}")

### Ceiling threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts ceiling threshold: {adata.shape}")
    total_counts_ceiling = adata.obs['total_counts'].quantile(0.95)
    adata_list[i] = adata[adata.obs['total_counts'] < total_counts_ceiling, :]
    print(f"Dataset {i} shape after total_counts ceiling threshold: {adata_list[i].shape}")

## Number of genes per cell
### Floor threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before n_genes_by_counts floor threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] > 500, :].copy()
    print(f"Dataset {i} shape after n_genes_by_counts floor threshold: {adata_list[i].shape}")

### Ceiling threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before n_genes_by_counts ceiling threshold: {adata.shape}")
    ngenes_ceiling = adata.obs['n_genes_by_counts'].quantile(0.97)
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] < ngenes_ceiling, :]
    print(f"Dataset {i} shape after n_genes_by_counts ceiling threshold: {adata_list[i].shape}")
    
## Ribosomal counts
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before pct_counts_ribosomal ceiling threshold: {adata.shape}")
    ribo_ceiling = adata.obs['pct_counts_ribosomal'].quantile(0.99)
    adata_list[i] = adata[adata.obs['pct_counts_ribosomal'] < ribo_ceiling, :]
    print(f"Dataset {i} shape after pct_counts_ribosomal ceiling threshold: {adata_list[i].shape}")

sample_A, sample_B, sample_C, sample_D, sample_E, sample_F = adata_list

## Mitochondrial percentage filter
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before mitochondrial filter threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs.pct_counts_mt < 1, :].copy()
    print(f"Dataset {i} shape after mitochondrial filter threshold: {adata.shape}")

sample_A, sample_B, sample_C, sample_D, sample_E, sample_F = adata_list

###############################################################################

# RECALCULATE QC METRICS

for i, adata in enumerate(adata_list):
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)

## Print QC summary data
qc_metrics = ['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal', 'n_genes_by_counts']

for i, adata in enumerate(adata_list, start=1):
    print(f"Dataset {i}:")
    for metric in qc_metrics:
        mean_value = adata.obs[metric].mean()
        print(f"Mean {metric}: {mean_value}")
    print("-" * 30) 
    
###############################################################################

# VISUALIZE EFFECT OF QC FILTERING

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal'])


os.chdir('/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/qc/filtered')

## Map each AnnData object to the specific group identifier
identifiers = ['A', 'B', 'C', 'D', 'E', 'F']
adata_list = [sample_A, sample_B, sample_C, sample_D, sample_E, sample_F]

for adata, id in zip(adata_list, identifiers):
    keys = ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal']
    for key in keys:
        sc.pl.violin(adata, key, show=False)
        plt.savefig(f"{key}_violin_group_{id}_filtered.pdf")
        plt.close()  
        
###############################################################################

# DOUBLET DETECTION

def run_scrublet(adata_list):

    scrublet_rows = []

    for i, adata in enumerate(adata_list):
        ## Set up Scrublet
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.10, random_state = 0)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                                  min_cells=3, 
                                                                  min_gene_variability_pctl=85, 
                                                                  n_prin_comps=50)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets

        ## Plot Histograms
        plt.figure(figsize=(10, 6))
        sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
        sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
        plt.title(f'Scrublet Doublet Score Distribution for Sample {i+1}')
        plt.xlabel('Doublet Score')
        plt.ylabel('Density')
        plt.legend()
        plt.grid(False)
        plt.show()

        ## Extract barcodes
        cell_barcodes = adata.obs.index

        ## Store Scrublet data with barcodes for each sample in list
        for barcode, obs_score, sim_score, pred_doublet in zip(cell_barcodes, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_, predicted_doublets):
            scrublet_rows.append({'Sample_Index': i+1, 
                                  'Cell_Barcode': barcode,
                                  'Observed_Score': obs_score, 
                                  'Simulated_Score': sim_score, 
                                  'Predicted_Doublet': pred_doublet})

    ## Create df from list of rows
    scrublet_df = pd.DataFrame(scrublet_rows)
    return scrublet_df

## Run Scrublet
scrublet_result = run_scrublet(adata_list)


###############################################################################

# VISUALIZE DOUBLET DETECTION RESULTS

## Define identifiers for the datasets
identifiers = ['A', 'B', 'C', 'D', 'E', 'F']

## Initialize empty list
data_for_plot = []

## Count singlets and doublets in each sample
for adata, id in zip(adata_list, identifiers):

    singlets_count = (adata.obs['predicted_doublets'] == False).sum()
    doublets_count = (adata.obs['predicted_doublets'] == True).sum()
    data_for_plot.append([id, singlets_count, doublets_count])


df = pd.DataFrame(data_for_plot, columns=['Dataset', 'Singlets', 'Doublets'])

## Plot doublet detection results
fig, ax = plt.subplots()
df.plot(x='Dataset', y=['Singlets', 'Doublets'], kind='bar', ax=ax, color=['blue', 'red'])
ax.set_title('Doublet detection results')
ax.set_xlabel('Dataset')
ax.set_ylabel('Count')
plt.xticks(rotation=0)
plt.show()

###############################################################################

# FILTER DOUBLETS

def filter_doublets(adata_list):
    filtered_list = []
    for i, adata in enumerate(adata_list):
        print(f"Dataset {i + 1} shape before doublet filter: {adata.shape}")
        filtered_adata = adata[adata.obs['predicted_doublets'] == False].copy()
        filtered_list.append(filtered_adata)
    return filtered_list

## Perform filtering
adata_list = filter_doublets(adata_list)

## Update list reference
sample_A, sample_B, sample_C, sample_D, sample_E, sample_F = adata_list

for i, adata in enumerate(adata_list):
    print(f"Dataset {i + 1} shape after doublet filter: {adata.shape}")


###############################################################################

# CONCATENATION

## Perform concatenation
adata_concat = sample_A.concatenate(sample_B, sample_C, sample_D, sample_E, sample_F, batch_key='group', batch_categories=['A', 'B', 'C', 'D', 'E', 'F'] )
adata_concat.obs

###############################################################################

# FINISH PRE-PROCESSING

## Filter genes appearing in few cells in concatenated dataset
sc.pp.filter_genes(adata_concat, min_cells=3)

## Normalization and Log1p transformation after concatenation
adata_concat.layers['counts'] = adata_concat.X.copy()
sc.pp.normalize_total(adata_concat) # median-based normalization
adata_concat.layers['normalized'] = adata_concat.X.copy()
sc.pp.log1p(adata_concat)
adata_concat.layers['log1p'] = adata_concat.X.copy()
adata_concat.raw = adata_concat

## Create separate object for full gene set
bdata = adata_concat.copy()

###############################################################################

# FEATURE SELECTION

## Select Highly Variable Genes (HVGs)
sc.pp.highly_variable_genes(adata_concat, n_top_genes = 3000, subset = True, layer = 'counts', flavor = "seurat_v3") # raw counts for seurat-v3

## Scale data
sc.pp.scale(adata_concat)
adata_concat.layers['scaled'] = adata_concat.X.copy()

###############################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata_concat.X = csr_matrix(adata_concat.X)
bdata.X = csr_matrix(bdata.X)

## save HVG version for scVI training
adata_concat.write_h5ad(os.path.join(save_dir, '1-tbi-seq-hvg.h5ad'))

## save full genome version for differential expression
bdata.write_h5ad(os.path.join(save_dir, '2-tbi-seq-full.h5ad'))

###############################################################################