## Background

### Experimental design
This advanced example workflow uses a snRNA-seq dataset of mouse hippocampus. PS19 mice (model for frontotemporal dementia) were subjected to mild traumatic brain injury (TBI), which is known to exacerbate  tau pathology in humans. An experimental AAV is used to rescue pathology and sequencing reveals changes to different cell populations, with an emphasis on neuronal populations and transcriptional signatures associated with pathology vs. health.

The workflow in this folder is adapted from an iterative workflow comprising of multiple phases of analysis, annotation, and visualization in preparation for manuscript submission.

---

## Scripts

### `01_pre-processing.py`
Performs quantile-based QC filtering to remove low-quality cells, detect doublets using Scrublet, and apply ceiling/floor thresholds for total and gene counts, ribosomal, and mitochondrial content. The cleaned data is normalized, highly variable genes (HVGs) are selected, and the processed dataset is exported for further analysis.

### `02_train-multi-models.ipynb`
Trains eight different scVI models based on varying filtering and doublet exclusion criteria, utilizing HVGs to control noise in UMAP visualization. After model selection and clustering, these models are integrated with full-genome data for downstream differential expression analysis.

### `03_eda-explore-model-results.py`
Conducts exploratory data analysis (EDA) on snRNA-seq data, including SCVI model visualization, UMAP plotting, clustering, and QC metrics evaluation. Outputs include CSV reports summarizing cluster statistics, QC data, and doublet scores across models.

### `04_eda-explore-data.py`
Explores gene expression patterns and QC metrics in the TBI dataset, generating UMAP visualizations for key genes and treatment groups, while analyzing doublet scores and cell distributions across treatment conditions.

### `05_annotate-clusters.py`
Processes TBI snRNA-seq data with SCVI to generate UMAPs, annotate clusters by cell types, and visualize gene markers. The script outputs summary plots and annotated datasets for further downstream analysis.

### `06_merge-hvg-full.py`
Integrates HVG annotations and full-genome data by matching cell barcodes between two datasets. Relevant metadata is transferred, and the updated full-genome dataset is saved for subsequent differential expression testing.

### `07_differential-expression.py`
Performs Wald test-based differential expression analysis across cell types and treatment groups, filtering results and generating heatmaps of up- and downregulated genes, enabling visual exploration of cell-type-specific gene expression changes.

### `08_data-visualization.py`
Generates key visualizations for TBI analysis figures, including UMAPs, differential expression heatmaps, cell type compositions, and gene expression plots related to neuronal health and pathology. Visualizations are tailored for publication-ready figures.

### `09_export-anndata-components.py`
Prepares snRNA-seq data for R/Seurat interoperability by cleaning annotations, recalibrating clusters, and exporting data layers (counts, log1p, normalized matrices) alongside dimensional reductions and metadata. Outputs ensure smooth cross-environment analysis.

### `10_seurat-interoperability.R`
Manually constructs Seurat objects from exported AnnData data, integrating raw, normalized, and log-transformed counts with metadata and dimensional reductions (PCA, UMAP, latent embeddings). Produces full-genome and HVG Seurat objects for visualization and downstream analysis in R.

---

## Data availability
This section will be updated soon with a link to the dataset in GEO.