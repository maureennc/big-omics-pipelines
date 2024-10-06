#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 14:52:28 2024

@author: maureen
"""

import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

################################################################################################################################

# IMPORT DATA
data_dir = "/Users/maureen/Documents/projects/harris-lab/Isaac/hunter_dataset/revisions/gene-ontology"

# GO BP
bp = pd.read_csv(os.path.join(data_dir, 'top_50_BP_DIRECT_pathways.csv'))
bp_sorted = bp.sort_values(by='Fold Enrichment', ascending=False)
bp_sorted['Term'] = bp_sorted['Term'].str.split(':').str[1]
bp_sorted['Term'] = bp_sorted['Term'].str.split('~').str[1]
bp_sorted['-log10(PValue)'] = -np.log10(bp_sorted['PValue'])

## KEGG

kegg = pd.read_csv(os.path.join(data_dir, 'top_50_kegg_pathways.csv'))
kegg_sorted = kegg.sort_values(by='Fold Enrichment', ascending=False)
kegg_sorted['Term'] = kegg_sorted['Term'].str.split(':').str[1]
kegg_sorted['-log10(PValue)'] = -np.log10(kegg_sorted['PValue'])

save_dir = "/Users/maureen/Documents/projects/harris-lab/Isaac/hunter_dataset/figures/2024-05-18"

################################################################################################################################

# GO BP ENRICHMENT

## Top 50
norm = mcolors.Normalize(vmin=bp_sorted['-log10(PValue)'].min(), vmax=bp_sorted['-log10(PValue)'].max())


fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
colors = plt.cm.viridis(norm(bp_sorted['-log10(PValue)'])).tolist()  # Convert to list
sns.barplot(x='Fold Enrichment', y='Term', data=bp_sorted, palette=colors, ax=ax)

ax.set_xlabel('Fold Enrichment', fontsize=12)
ax.set_ylabel('')
ax.set_title('Top 50 GO BP Terms Enriched in Casp1+Cx3cr1+ cells', fontsize=14)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

# Create a color bar
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, label='-log10(PValue)', aspect=20, shrink=0.5)
cbar.set_label('-log10(PValue)', size=12)
cbar.ax.tick_params(labelsize=10)
tick_values = np.linspace(start=bp_sorted['-log10(PValue)'].min(), stop=bp_sorted['-log10(PValue)'].max(), num=5)
cbar.set_ticks(tick_values)
cbar.set_ticklabels(['{:.1f}'.format(x) for x in -np.log10(np.power(10, -tick_values))])
plt.tight_layout()

plt.savefig(os.path.join(save_dir, 'go-top-50.pdf'), format='pdf', bbox_inches='tight')
plt.show()



## Top 25
top_25_bp = bp_sorted.head(25)
norm = mcolors.Normalize(vmin=top_25_bp['-log10(PValue)'].min(), vmax=top_25_bp['-log10(PValue)'].max())


fig, ax = plt.subplots(figsize=(5, 12), dpi=300)
colors = plt.cm.viridis(norm(top_25_bp['-log10(PValue)'])).tolist()  # Convert to list
barplot = sns.barplot(x='Fold Enrichment', y='Term', data=top_25_bp, palette=colors, ax=ax)

ax.set_xlabel('Fold Enrichment', fontsize=10)
ax.set_ylabel('')
ax.set_title('Top 25 GO BP Terms Enriched in Casp1+Cx3cr1+ cells', fontsize=14, pad = 20)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

# Wrap the y-axis labels to make the plot less wide
ax.set_yticklabels(ax.get_yticklabels(), wrap=True, fontsize = 10)

# Create a color bar
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, label='-log10(PValue)', aspect=10, shrink=0.5)
cbar.set_label('-log10(PValue)', size=12)
cbar.ax.tick_params(labelsize=10)
tick_values = np.linspace(start=top_25_bp['-log10(PValue)'].min(), stop=top_25_bp['-log10(PValue)'].max(), num=5)
cbar.set_ticks(tick_values)
cbar.set_ticklabels(['{:.1f}'.format(x) for x in -np.log10(np.power(10, -tick_values))])

plt.tight_layout()

plt.savefig(os.path.join(save_dir, 'go-top-25.pdf'), format='pdf', bbox_inches='tight')
plt.show()


## Top 15
top_15_bp = bp_sorted.head(15)
norm = mcolors.Normalize(vmin=top_15_bp['-log10(PValue)'].min(), vmax=top_15_bp['-log10(PValue)'].max())


fig, ax = plt.subplots(figsize=(5, 12))
colors = plt.cm.viridis(norm(top_15_bp['-log10(PValue)'])).tolist()  # Convert to list
barplot = sns.barplot(x='Fold Enrichment', y='Term', data=top_15_bp, palette=colors, ax=ax)

ax.set_xlabel('Fold Enrichment', fontsize=14)
ax.set_ylabel('')
ax.set_title('Top 15 GO BP Terms Enriched in Casp1+Cx3cr1+ cells', fontsize=14, pad = 20)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

# Wrap the y-axis labels to make the plot less wide
ax.set_yticklabels(ax.get_yticklabels(), wrap=True, fontsize = 12)

# Create a color bar
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, label='-log10(PValue)', aspect=8, shrink=0.5)
cbar.set_label('-log10(PValue)', size=12)
cbar.ax.tick_params(labelsize=10)
tick_values = np.linspace(start=top_15_bp['-log10(PValue)'].min(), stop=top_15_bp['-log10(PValue)'].max(), num=5)
cbar.set_ticks(tick_values)
cbar.set_ticklabels(['{:.1f}'.format(x) for x in -np.log10(np.power(10, -tick_values))])

plt.tight_layout()
plt.grid(False)

plt.savefig(os.path.join(save_dir, 'go-top-15.pdf'), format='pdf', bbox_inches='tight')
plt.show()

################################################################################################################################

# KEGG ENRICHMENT

## Top 50
norm = mcolors.Normalize(vmin=kegg_sorted['-log10(PValue)'].min(), vmax=kegg_sorted['-log10(PValue)'].max())

fig, ax = plt.subplots(figsize=(5, 12), dpi=300)
colors = plt.cm.viridis(norm(kegg_sorted['-log10(PValue)'])).tolist()
sns.barplot(x='Fold Enrichment', y='Term', data=bp_sorted, palette=colors, ax=ax)

ax.set_xlabel('Fold Enrichment', fontsize=10)
ax.set_ylabel('')
ax.set_title('Top 50 KEGG Pathways Enriched in Casp1+Cx3cr1+ cells', fontsize=14)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

# Create a color bar
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, label='-log10(PValue)', aspect=8, shrink=0.5)
cbar.set_label('-log10(PValue)', size=12)
cbar.ax.tick_params(labelsize=10)
tick_values = np.linspace(start=kegg_sorted['-log10(PValue)'].min(), stop=kegg_sorted['-log10(PValue)'].max(), num=5)
cbar.set_ticks(tick_values)
cbar.set_ticklabels(['{:.1f}'.format(x) for x in -np.log10(np.power(10, -tick_values))])
plt.tight_layout()
plt.grid(False)

plt.savefig(os.path.join(save_dir, 'kegg-top-50.pdf'), format='pdf', bbox_inches='tight')
plt.show()



## Top 25
top_25_kegg = kegg_sorted.head(25)
norm = mcolors.Normalize(vmin=top_25_kegg['-log10(PValue)'].min(), vmax=top_25_kegg['-log10(PValue)'].max())


fig, ax = plt.subplots(figsize=(5, 12), dpi=300)
colors = plt.cm.viridis(norm(top_25_kegg['-log10(PValue)'])).tolist()  # Convert to list
barplot = sns.barplot(x='Fold Enrichment', y='Term', data=top_25_kegg, palette=colors, ax=ax)

ax.set_xlabel('Fold Enrichment', fontsize=10)
ax.set_ylabel('')
ax.set_title('Top 25 KEGG Terms Enriched in Casp1+Cx3cr1+ cells', fontsize=14, pad = 20)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

# Wrap the y-axis labels to make the plot less wide
ax.set_yticklabels(ax.get_yticklabels(), wrap=True, fontsize = 10)

# Create a color bar
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, label='-log10(PValue)', aspect=8, shrink=0.5)
cbar.set_label('-log10(PValue)', size=12)
cbar.ax.tick_params(labelsize=10)
tick_values = np.linspace(start=top_25_kegg['-log10(PValue)'].min(), stop=top_25_kegg['-log10(PValue)'].max(), num=5)
cbar.set_ticks(tick_values)
cbar.set_ticklabels(['{:.1f}'.format(x) for x in -np.log10(np.power(10, -tick_values))])

plt.tight_layout()
plt.grid(False)

plt.savefig(os.path.join(save_dir, 'kegg-top-25.pdf'), format='pdf', bbox_inches='tight')
plt.show()



## Top 15
top_15_kegg = kegg_sorted.head(15)
norm = mcolors.Normalize(vmin=top_15_kegg['-log10(PValue)'].min(), vmax=top_15_kegg['-log10(PValue)'].max())


fig, ax = plt.subplots(figsize=(5, 12))
colors = plt.cm.viridis(norm(top_15_kegg['-log10(PValue)'])).tolist()  # Convert to list
barplot = sns.barplot(x='Fold Enrichment', y='Term', data=top_15_kegg, palette=colors, ax=ax)

ax.set_xlabel('Fold Enrichment', fontsize=14)
ax.set_ylabel('')
ax.set_title('Top 15 KEGG Pathways Enriched in Casp1+Cx3cr1+ cells', fontsize=14, pad = 20)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

# Wrap the y-axis labels to make the plot less wide
ax.set_yticklabels(ax.get_yticklabels(), wrap=True, fontsize = 12)

# Create a color bar
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, label='-log10(PValue)', aspect=8, shrink=0.5)
cbar.set_label('-log10(PValue)', size=12)
cbar.ax.tick_params(labelsize=10)
tick_values = np.linspace(start=top_15_kegg['-log10(PValue)'].min(), stop=top_15_kegg['-log10(PValue)'].max(), num=5)
cbar.set_ticks(tick_values)
cbar.set_ticklabels(['{:.1f}'.format(x) for x in -np.log10(np.power(10, -tick_values))])

plt.tight_layout()
plt.grid(False)

plt.savefig(os.path.join(save_dir, 'kegg-top-15.pdf'), format='pdf', bbox_inches='tight')
plt.show()


################################################################################################################################