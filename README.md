# About

This code repository has been created as a resource at the Center for Brain Immunology and Glia (BIG) and Neuroscience department at the University of Virginia. BIG is a community of biomedical researchers dedicated to the investigation of inflammatory processes across neurological conditions including infection and neurodegeneration. This repository houses a collection of data pipelines and example workflows put together to support the  diverse transcriptomic and spatial biology experiments performed in the BIG Center.

![MERFISH UMAP](visualization/figures/merfish-spatial-scatter.png)

# Data processing pipelines
## 1. `bulk-RNA-seq-data`
Start with .fastq files from paired-end sequencing and run a trimmomatic-salmon pipeline. Read data into an R environment and proceed with a DESeq2-driven analysis. Also includes code for gene ontology overrepresentation analysis for differential expression results. Featured dataset is a sequencing experiment from control and T. gondii-infected brains. The featured dataset (Harris lab) was generated to obtain infection-specific FPKM (abundance) values to guide in the creation of 500 and 1000-plex MERFISH panels, with the goal of preventing optical crowding during data generation. 

## 2. `merfish-data`
Perform segmentation, data processing, and computational analysis on MERFISH data collected from control and T. gondii-infected mouse brains (Harris lab). Segmentation is performed on the Rivanna/Afton HPC using the cellpose 2.0 cyto2 algorithm via the Vizgen post-processing tool (VPT). After segmentation, transcripts are partitioned into cell boundaries. The dataset is imported into a Python environment and assembled into an AnnData object for single-cell analysis. See [poster PDF](visualization/figures/MERFISH_HPC_Pipeline_Cowan_RCSymposium2024_poster.pdf) for a comprehensive overview of the computational workflow.
   
## 3. `nanostring-cosmx-data`
Prepare and analyze Nanostring CosMx SMI data. Example workflow features Nanostring demo data and a mouse brain dataset from an aging-associated Neuro-COVID19 project (Lukens lab). CosMx data is pre-processed using AtoMx software wth cellpose segmentation prior to transfer to an AWS S3 bucket for subsequent processing using a lab-specific cloud-computing infrastructure. Scripting for data processing 
   
## 4. `nanostring-geomx-data`
Prepare and analyze Nanostring GeoMx Digital Satial Profiler (DSP) data. GeoMx is an ROI-based spatial data is analyzed using the `GeoMx tools` bioconductor package. See [vignette](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html#4_QC__Pre-processing). Featured dataset generated from the mouse olfactory system during SARS-CoV-2 infection (Lukens lab) and analysis performed in R.
   
## 5. `scRNA-seq-data`
Run cellranger and perform single-cell analysis on 10x genomics data. Dataset features immune cells that were FACS-sorted from T. gondii-infected mouse brains (Harris lab). Data cleaning involves filtering on QC parameters using a dynamic quantile approach that scales to each sample, and scrublet for doublet detection. This section includes example scripts for cell type annotation and identification and differential expression. Trajectory inference analysis (RNA velocity) is performed to examine the microglial transition from homeostatic to a neurodegeneration-associated transcriptional state during parasitic infection using tools including samtools, velocyto and scVelo.  
   
## 6. `snRNA-seq-data`
Includes basic and in-depth analysis workflows using two datasets using genetic mouse models of Alzheimer's Disease (Lukens lab). The workflow for single nuclei RNA-sequencing data is very similar to single-cell, with additional considerations for QC parameters such as increased sparsity and lower mitochondrial read fraction.
   

# Other sections
## 1. `tools`
This section houses scripts for a collection of bioinformatics tools for single-cell analysis,  including cellranger, velocyto, and scvi-tools.
   
## 2. `visualization`: 
This section is subdivided into scripts, instructions, and figures for high-dimensional data visualization.



