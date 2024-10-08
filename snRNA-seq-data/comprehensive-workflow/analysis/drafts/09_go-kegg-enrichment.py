#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 17:27:55 2024

@author: maureen

GO Enrichment plotting

"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap

###############################################################################

# GENERATE GENE LISTS FOR GO OVERREPRESENTATION ANALYSIS

de_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de/cell_type/filtered'

ce_results = pd.read_csv(os.path.join(de_dir, 'CE-comparison-cell_type-filtered.csv'))

## Create subsets
ce_excitatory = ce_results[ce_results['cell_type'] == 'Excitatory neuron'].copy()
ce_inhibitory = ce_results[ce_results['cell_type'] == 'Inhibitory neuron'].copy()

## Set cutoffs
ce_excitatory = ce_excitatory[(abs(ce_excitatory['log2fc']) > 0.1) & (ce_excitatory['qval'] <= 0.05)].copy()
ce_inhibitory = ce_inhibitory[(abs(ce_inhibitory['log2fc']) > 0.1) & (ce_inhibitory['qval'] <= 0.05)].copy()

## Sort by log2fc and adjusted p value
ce_excitatory_sorted = ce_excitatory.sort_values(by=['log2fc', 'qval'], ascending=[False, True])
ce_inhibitory_sorted = ce_inhibitory.sort_values(by=['log2fc', 'qval'], ascending=[False, True])

###############################################################################

# CREATE GENE LISTS

input_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/gsea/CE/input'

ce_excitatory_sorted.to_csv(os.path.join(input_dir, 'CE-excitatory-neuron-degs-list.csv'))
ce_inhibitory_sorted.to_csv(os.path.join(input_dir, 'CE-inhibitory-neuron-degs-list.csv'))

#v1 is with a log2fc of 0.1 -- used for GSEA
#v2 is with a log2fc of 0.5 -- insufficient genes

###############################################################################

# VISUALIZE RESULTS

output_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/gsea/CE/output'

## Import results
excitatory_results = pd.read_csv(os.path.join(output_dir, 'GO-BP-excitatory-group-E-enriched.txt'), delimiter = '\t')
inhibitory_results = pd.read_csv(os.path.join(output_dir, 'GO-BP-inhibitory-group-E-enriched.txt'), delimiter = '\t')

## Create -log10pvalue column
excitatory_results['-log10 pval'] = -np.log10(excitatory_results['Adjusted P-value'])
inhibitory_results['-log10 pval'] = -np.log10(inhibitory_results['Adjusted P-value'])


## Convert overlap column
inhibitory_results[['Overlap_count', 'GO_term_total']] = inhibitory_results['Overlap'].str.split('/', expand=True).astype(int)
inhibitory_results['Gene Ratio'] = inhibitory_results['Overlap_count'] / inhibitory_results['GO_term_total']

excitatory_results[['Overlap_count', 'GO_term_total']] = excitatory_results['Overlap'].str.split('/', expand=True).astype(int)
excitatory_results['Gene Ratio'] = excitatory_results['Overlap_count'] / excitatory_results['GO_term_total']

excitatory_terms = excitatory_results.head(10)
inhibitory_terms = inhibitory_results.head(10)

###############################################################################

# PLOTTING 

gsea_save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/go/CE'


def wrap_labels(labels, width):
    return [textwrap.fill(label, width) for label in labels]

excitatory_terms['Wrapped Term'] = wrap_labels(excitatory_terms['Term'], width=40)
inhibitory_terms['Wrapped Term'] = wrap_labels(inhibitory_terms['Term'], width=40)


## Plotting excitatory neurons
plt.figure(figsize=(10, 9))
sns.barplot(
    data=excitatory_terms,
    x='-log10 pval',
    y='Wrapped Term',
    color='skyblue'
)
plt.xlabel('-log10 pval (adjusted)')
plt.ylabel('')
plt.title('')
plt.tight_layout()
plt.savefig(os.path.join(gsea_save_dir, 'excitatory-neurons-CE-GO-top10.pdf'))
plt.show()


## Plotting inhibitory neurons
plt.figure(figsize=(10, 9))
sns.barplot(
    data=inhibitory_terms,
    x='-log10 pval',
    y='Wrapped Term',
    color='skyblue'
)
plt.xlabel('-log10 pval (adjusted)')
plt.ylabel('')
plt.title('')
plt.tight_layout()
plt.savefig(os.path.join(gsea_save_dir, 'inhibitory-neurons-CE-GO-top10.pdf'))
plt.show()

###############################################################################