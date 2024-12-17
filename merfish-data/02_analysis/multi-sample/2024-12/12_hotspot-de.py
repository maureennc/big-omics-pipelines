#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:26:16 2024

@author: maureen
"""

###############################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/processed'

adata = sc.read_h5ad(os.path.join(data_dir, 'concat-split-trained-annotated.h5ad'))

## Prepare data

adata.X = adata.layers['log1p'].copy()
### Infected only
adata = adata[adata.obs['condition'] == 'Infected'].copy()

### convert from boolean to string
adata.obs['inflammatory_foci'] = adata.obs['inflammatory_foci'].astype(str)

print(adata.obs['inflammatory_foci'].unique()) 


###############################################################################

# DE

adata.obs['inflammatory_foci'] = adata.obs['inflammatory_foci'].astype(str)

# Function to compute DE for inflammatory foci vs. non-foci
def compute_de_foci(adata, min_cells=3):
    results = []  # Store DE results

    for cell_type in adata.obs['cell_type'].unique():
        adata_subset = adata[adata.obs['cell_type'] == cell_type].copy()
        sc.pp.filter_genes(adata_subset, min_cells=min_cells)  # Filter genes

        group_counts = adata_subset.obs['inflammatory_foci'].value_counts()
        if 'True' in group_counts and 'False' in group_counts and \
           group_counts['True'] > 40 and group_counts['False'] > 40:

            print(f"Comparing foci vs non-foci for: {cell_type}")
            sc.tl.rank_genes_groups(
                adata_subset, groupby='inflammatory_foci', reference='False',
                method='wilcoxon', pts=True
            )

            # Extract DE results
            de_results = sc.get.rank_genes_groups_df(adata_subset, group='True')
            comparison_nonzero = (adata_subset[adata_subset.obs['inflammatory_foci'] == 'True'].X > 0).mean(axis=0).A1
            reference_nonzero = (adata_subset[adata_subset.obs['inflammatory_foci'] == 'False'].X > 0).mean(axis=0).A1

            de_results['pct_nz_group'] = comparison_nonzero
            de_results['pct_nz_reference'] = reference_nonzero
            de_results['pct_nz_mean'] = (comparison_nonzero + reference_nonzero) / 2
            de_results['cell_type'] = cell_type

            results.append(de_results)

    return pd.concat(results, ignore_index=True)

# Compute DE analysis
foci_vs_nonfoci_df = compute_de_foci(adata)


###############################################################################

# SUMMARY HEATMAP

logfc_threshold = 0.5
qval_threshold = 0.01
mean_expression_threshold = 0.10

# Filter results based on thresholds
filtered_results = foci_vs_nonfoci_df[
    (abs(foci_vs_nonfoci_df['logfoldchanges']) > logfc_threshold) &
    (foci_vs_nonfoci_df['pvals_adj'] < qval_threshold) &
    (foci_vs_nonfoci_df['pct_nz_mean'] > mean_expression_threshold)
]

# Initialize DataFrames for DEG counts
order = [
    'Astrocyte', 'Choroid plexus', 
    'Excitatory neuron', 'Inhibitory neuron',
    'OPC', 'Oligodendrocyte', 'CD4+ T cell', 'CD8+ T cell', 'Macrophage',  'Microglia'
] # Remove vascular types 'Endothelial cell' 'Pericyte'



cell_type_reference = pd.DataFrame(order, columns=['cell_type'])
deg_counts_all = cell_type_reference.copy()
deg_counts_all['Upregulated'] = 0
deg_counts_all['Downregulated'] = 0

# Count DEGs for each cell type
for cell_type in filtered_results['cell_type'].unique():
    df_cell_type = filtered_results[filtered_results['cell_type'] == cell_type]
    upregulated = df_cell_type[df_cell_type['logfoldchanges'] > 0].shape[0]
    downregulated = df_cell_type[df_cell_type['logfoldchanges'] < 0].shape[0]

    deg_counts_all.loc[deg_counts_all['cell_type'] == cell_type, 'Upregulated'] = upregulated
    deg_counts_all.loc[deg_counts_all['cell_type'] == cell_type, 'Downregulated'] = downregulated

# Prepare DataFrame for heatmap
deg_counts_all = deg_counts_all[['cell_type', 'Upregulated', 'Downregulated']]
deg_counts_all = deg_counts_all.set_index('cell_type').T  # Transpose
deg_counts_all.index = ['UP', 'DOWN']

# Plot horizontal heatmap
plt.figure(figsize=(5, 2.5))
sns.heatmap(deg_counts_all, annot=True, cmap='rocket', fmt='g', cbar_kws={'label': '# DEGs'})
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

###############################################################################

# EXPORT

save_dir = '/Users/maureen/Desktop/toxo-merfish/spreadsheets'

naive_vs_infected_df.to_csv(os.path.join(save_dir, 'naive-vs-infected-de.csv'), index = False)
foci_vs_nonfoci_df.to_csv(os.path.join(save_dir, 'foci-vs-non-foci-de.csv'), index = False)

###############################################################################