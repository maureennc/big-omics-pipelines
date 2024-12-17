#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:00:06 2024

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
import squidpy as sq
import random

###############################################################################

# SETTINGS


## Matplotlib
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 500

###############################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/processed'

adata = sc.read_h5ad(os.path.join(data_dir, 'concat-split-trained-annotated.h5ad'))

###############################################################################

# PIE CHARTS

## Hematopoietic vs. non-hematopoietic
hematopoietic = ["Microglia", "CD4+ T cell", "CD8+ T cell"]
non_hematopoietic = [
    "Astrocyte", "Excitatory neuron", "Inhibitory neuron", 
    "OPC", "Oligodendrocyte", "Endothelial cell", 
    "Pericyte", "Choroid plexus"
]

# Add a new column in adata.obs to classify each cell
adata.obs['compartment'] = adata.obs['cell_type'].apply(
    lambda x: 'Hematopoietic' if x in hematopoietic else 'Non-Hematopoietic'
)

# Calculate the fraction of each compartment within Naive and Infected groups
compartment_counts = (
    adata.obs.groupby(['condition', 'compartment']).size().reset_index(name='count')
)

# Calculate the total number of cells per condition
total_counts = adata.obs.groupby('condition').size().reset_index(name='total')

# Merge the total counts with the compartment counts
compartment_counts = pd.merge(compartment_counts, total_counts, on='condition')

# Calculate the fraction of each compartment per condition
compartment_counts['fraction'] = compartment_counts['count'] / compartment_counts['total']

# Inspect the resulting DataFrame
print(compartment_counts)


## Plot
import matplotlib.pyplot as plt
import seaborn as sns

# Get two colors from the 'rocket' colormap
rocket_colors = sns.color_palette('rocket', 2)
compartment_colors = {'Hematopoietic': rocket_colors[0], 'Non-Hematopoietic': rocket_colors[1]}

# Function to format the percentage labels
def func(pct, threshold=5):
    """Format the label to only show percentage if above the threshold."""
    if pct > threshold:
        return f"{pct:.1f}%"  # Display percentage only
    else:
        return ""  # No label if the slice is <= threshold

# Function to plot a pie chart for a specific condition
def plot_pie(condition):
    data = compartment_counts[compartment_counts['condition'] == condition]
    
    fig, ax = plt.subplots(figsize=(5, 5))  # Smaller pie chart size

    # Create the pie chart with conditional percentage labels
    wedges, texts, autotexts = ax.pie(
        data['fraction'], 
        startangle=90,
        colors=[compartment_colors[c] for c in data['compartment']],
        explode=(0.1,) * len(data),  # Slightly explode slices
        pctdistance=0.6,  # Adjust position of percentage labels
        autopct=lambda pct: func(pct, threshold=5)  # Apply threshold logic
    )
    
    # Customize the text appearance
    for autotext in autotexts:
        autotext.set_fontsize(16)  # Increase font size for percentages
        autotext.set_color('black')  # Set color to black

    ax.axis('equal')  # Ensure the pie chart is circular

    # Create a legend to the right of the pie chart
    handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=12) 
        for color in rocket_colors
    ]
    ax.legend(
        handles, ['Hematopoietic', 'Non-Hematopoietic'], 
        title="Compartments", fontsize=14, title_fontsize=16, loc="center left", 
        bbox_to_anchor=(1, 0.5)  # Move legend to the right
    )

    plt.show()

# Plot the two pie charts separately
plot_pie('Naive')
plot_pie('Infected')
