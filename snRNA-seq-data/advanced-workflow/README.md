---

# Advanced snRNA-seq Analysis Workflow

## Background

### Experimental Design
This advanced workflow focuses on a snRNA-seq dataset from the mouse hippocampus. PS19 mice (a model for frontotemporal dementia) were subjected to mild traumatic brain injury (TBI), known to exacerbate tau pathology in humans. An experimental AAV therapy was used to rescue pathology, and sequencing reveals transcriptional changes across cell populations, particularly in neuronal populations, highlighting signatures associated with disease vs. health.

The workflow is designed as an iterative process, involving multiple phases of analysis, annotation, and visualization to prepare for manuscript submission.

---

## Scripts

### `01_cellranger.sh`
Runs the Cell Ranger pipeline to process single-cell RNA-seq data. It specifies memory, cores, read lengths, and the reference transcriptome, automating the counting step for a specific sample to prepare it for downstream analysis.



### `02_pre-processing.py`
Applies quantile-based quality control (QC) filtering to remove low-quality cells, detect doublets using Scrublet, and apply thresholds for total/gene counts, ribosomal content, and mitochondrial content. The cleaned data is normalized, HVGs are selected, and the processed dataset is exported.



### `03_train-multi-models.ipynb`
Trains eight different scVI models based on varying filtering and doublet exclusion criteria, using HVGs to control noise in UMAP visualization. Models are then integrated with full-genome data for downstream differential expression analysis.



### `04_eda-explore-model-results.py`
Performs exploratory data analysis (EDA), visualizing SCVI models, UMAP plots, clustering, and QC metrics. Outputs include CSV reports summarizing cluster statistics, QC data, and doublet scores across multiple models.



### `05_eda-explore-data.py`
Analyzes gene expression patterns and QC metrics in the TBI dataset, generating UMAP visualizations for key genes and treatment groups, while assessing doublet scores and cell distribution across conditions.



### `06_annotate-clusters.py`
Processes TBI snRNA-seq data with SCVI to generate UMAPs, annotate clusters with cell types, and visualize gene markers. Annotated datasets and summary plots are exported for downstream analysis.



### `07_merge-hvg-full.py`
Integrates HVG annotations with full-genome data by matching cell barcodes between datasets. Metadata is transferred, and the updated full-genome dataset is saved for differential expression testing.



### `08_differential-expression.py`
Performs Wald test-based differential expression analysis across cell types and treatment groups, filtering results and generating heatmaps for visual exploration of gene expression changes by cell type.



### `09_data-visualization.py`
Generates key visualizations for the TBI analysis, including UMAPs, differential expression heatmaps, cell type compositions, and gene expression plots related to neuronal health and pathology. Visuals are prepared for publication-ready figures.



### `10_export-anndata-components.py`
Prepares snRNA-seq data for Seurat by exporting data layers (counts, log1p, normalized matrices) along with dimensional reductions and metadata. Ensures smooth interoperability between Python and R.



### `11_seurat-interoperability.R`
Constructs Seurat objects from the exported AnnData data, integrating raw, normalized, and log-transformed counts with metadata and dimensional reductions (PCA, UMAP, latent embeddings). Full-genome and HVG Seurat objects are created for visualization and downstream analysis in R.

---

## Data Availability
A link to the dataset in GEO will be added soon.

---
