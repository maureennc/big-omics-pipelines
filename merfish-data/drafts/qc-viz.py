#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 13:04:03 2024

@author: maureen
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MERFISH QC-viz
Created on Sun Dec 31 11:30:59 2023

@author: maureen
"""

# Starting point - Processing anndata lists
## import data

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

#adata_list = [adata1, adata2, adata3, adata4, adata5]

######################################################################

#QC1: PRE-FILTER QC
## apply to raw data

def plot_qc1(adata_list):

    qc_df = pd.DataFrame()  # Initialize the master DataFrame

    for i, adata in enumerate(adata_list):
        # QC Calculations
        sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)
        
        # Selecting relevant QC metrics
        qc_metrics = ['total_counts', 'n_genes_by_counts', 'volume']
        df_qc = adata.obs[qc_metrics].copy()
        df_qc = df_qc.melt(var_name='QC Metric', value_name='Value')
        df_qc['Sample'] = f'Sample {i+1}'

        # Append to the master DataFrame
        qc_df = pd.concat([qc_df, df_qc], ignore_index=True)

        # Plotting for the current sample
        g = sns.FacetGrid(df_qc, col='QC Metric', sharex=False, sharey=False)
        g.map_dataframe(sns.histplot, x='Value')

        # Adjusting subplot titles
        for ax in g.axes.flatten():
            ax.set_title(ax.get_title(), fontsize=10)  # Set smaller font size here

        g.fig.suptitle(f'QC for Sample {i+1}', fontsize=16)
        g.fig.subplots_adjust(top=0.85)
        plt.show()

    return qc_df

adata_list = [adata1, adata2, adata3, adata4, adata5
              ]
qc_df = plot_qc1(adata_list)
del qc_df #clear qc_df from memory


######################################################################

# PROCEED WITH FILTERING POOR QUALITY CELLS


######################################################################


# QC2: RUN TO PLOT INITIAL QC DATA AFTER FILTERING

def plot_qc2(adata_list):
    qc_df = pd.DataFrame()  # Initialize the master DataFrame

    for i, adata in enumerate(adata_list):
        # Recalculate QC Metrics after filtering
        sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)
        
        # Select relevant QC metrics
        qc_metrics = ['total_counts', 'n_genes_by_counts', 'volume']
        df_qc = adata.obs[qc_metrics].copy()
        df_qc = df_qc.melt(var_name='QC Metric', value_name='Value')
        df_qc['Sample'] = f'Sample {i+1}'

        # Append to master DataFrame
        qc_df = pd.concat([qc_df, df_qc], ignore_index=True)

        # Plotting for the current sample
        g = sns.FacetGrid(df_qc, col='QC Metric', sharex=False, sharey=False)
        g.map_dataframe(sns.histplot, x='Value')

        # Adjusting subplot titles
        for ax in g.axes.flatten():
            ax.set_title(ax.get_title(), fontsize=10)  # Set smaller font size here

        g.fig.suptitle(f'QC for Sample {i+1}', fontsize=16)
        g.fig.subplots_adjust(top=0.85)
        plt.show()

    return qc_df

qc_df = plot_qc2(adata_list)
del qc_df #clear qc_df from memory


################################################################


# QC3: CUSTOMIZED visualization on filtered cells


## QC for cells and genes

##################

### Cells


def create_obs_df(adata_list):
    obs_df = pd.DataFrame()
    
    for i, adata in enumerate(adata_list):
        # Extract metrics from .obs
        df_obs = adata.obs.copy()
        df_obs = df_obs.reset_index().melt(id_vars='index', var_name='Metric', value_name='Value')
        df_obs['Sample'] = f'Sample {i+1}'  # Corrected line with the closing quote
    
        # Append to the master DataFrame
        obs_df = pd.concat([obs_df, df_obs], ignore_index=True)
    
    return obs_df

obs_df = create_obs_df(adata_list)  # Correctly calling the function


##################

## Genes


def create_var_df(adata_list):
    var_df = pd.DataFrame()
    
    for i, adata in enumerate(adata_list):
        # Assuming 'mean_counts', 'pct_dropout_by_counts', 'n_cells_by_counts' are columns in adata.var
        df_var = adata.var[['mean_counts', 'pct_dropout_by_counts', 'n_cells_by_counts']].copy()
        df_var = df_var.reset_index().melt(id_vars='index', var_name='Metric', value_name='Value')
        df_var['Sample'] = f'Sample {i+1}'

        # Append to the master DataFrame
        var_df = pd.concat([var_df, df_var], ignore_index=True)
    
    return var_df

var_df = create_var_df(adata_list)

###################################


## QC3 Visualization

##################

## QC for Cells

def plot_obs_metrics(obs_df, adata_list):
    # Combine all metrics into one DataFrame for FacetGrid
    combined_df = pd.DataFrame()
    for i, _ in enumerate(adata_list):
        sample_name = f'Sample {i+1}'
        obs_sample_df = obs_df[obs_df['Sample'] == sample_name]
        combined_df = pd.concat([combined_df, obs_sample_df], ignore_index=True)

    # Define the metrics for plotting
    obs_histogram_metrics = ['transcript_count', 'total_counts', 'log1p_total_counts', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'volume', 'solidity', 'anisotropy', 'pct_counts_in_top_50_genes']

    # Create a FacetGrid for all metrics with three columns
    g = sns.FacetGrid(combined_df, col='Metric', hue='Sample', col_order=obs_histogram_metrics,
                      sharex=False, sharey=False, height=4, aspect=1.5, col_wrap=3)  # Set col_wrap to 3 for three columns
    g.map_dataframe(sns.histplot, 'Value', kde=True).add_legend()

    # Remove grid lines and set titles
    for ax in g.axes.flatten():
        ax.grid(False)  # Remove grid lines
        ax.set_title(ax.get_title().replace('Metric = ', ''), fontsize=10)

    plt.show()

# Example usage
plot_obs_metrics(obs_df, adata_list)




##################

## QC for Genes

def plot_var_metrics(var_df, adata_list):
    # Combine all metrics into one DataFrame for FacetGrid
    combined_df = pd.DataFrame()
    for i, _ in enumerate(adata_list):
        sample_name = f'Sample {i+1}'
        var_sample_df = var_df[var_df['Sample'] == sample_name]
        combined_df = pd.concat([combined_df, var_sample_df], ignore_index=True)

    # Define the metrics for plotting
    var_histogram_metrics = ['mean_counts', 'pct_dropout_by_counts', 'n_cells_by_counts']

    # Create a FacetGrid for all metrics with three columns
    g = sns.FacetGrid(combined_df, col='Metric', hue='Sample', col_order=var_histogram_metrics,
                      sharex=False, sharey=False, height=4, aspect=1.5, col_wrap=1)  # Set col_wrap to 3 for three columns
    g.map_dataframe(sns.histplot, 'Value', kde=True).add_legend()

    # Remove grid lines and set titles
    for ax in g.axes.flatten():
        ax.grid(False)  # Remove grid lines
        ax.set_title(ax.get_title().replace('Metric = ', ''), fontsize=10)

    plt.show()

# Example usage
plot_var_metrics(var_df, adata_list)



###################################################################

# RUN SCRUBLET


## Run scrublet

## Filter out cells by doublet threshold

####################################################################

# QC4: POST-DOUBLET REMOVAL
## examine effects of doublet removal


## UPDATE QC DATAFRAMES, obs_df and var_df

obs_df = create_obs_df(adata_list)  # Recompute obs metrics
var_df = create_var_df(adata_list)  # Recompute var metrics


##################

## VISUALIZE EFFECTS OF DOUBLET REMOVAL ON QC PLOTS FOR QC4

### Facet Grid (option 1)
plot_obs_metrics(obs_df, adata_list)
plot_var_metrics(var_df, adata_list)




##################

## PLOT NUMBER OF GENES PER COUNTS

def plot_n_genes_per_counts(adata_list):
    for i, adata in enumerate(adata_list):
        # Extract the 'n_genes_per_counts' data for the current AnnData object
        n_genes_data = adata.obs['n_genes_by_counts']

        # Create the plot
        plt.figure(figsize=(6, 4))
        sns.histplot(n_genes_data, kde=True)
        plt.title(f'Sample {i+1} - Histogram of n_genes_by_counts')
        plt.xlabel('n_genes_by_counts')
        plt.ylabel('Count')
        plt.show()

# Example usage
plot_n_genes_per_counts(adata_list)


######################################################################
