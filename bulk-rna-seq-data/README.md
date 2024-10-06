# Bulk-RNA-Sequencing
Bulk RNA Sequencing Pipeline for Harris lab naive vs. T. gondii-infected samples. Rivanna HPC used for bash and R scripting.

![Figure 2](../visualization/figures/bulk-rna-seq.png)

## Computational Workflow

### Step 1: Run `01_trimmomatic.sh`
Use Trimmomatic to trim out adapters from fastq files. Nextera adapters were used in example experiment. If unsure of which adapters should be trimmed, can run FastQC and read report.

Customizeable parameters used:
- Remove leading low quality or N bases (below quality 3) (LEADING:3)
- Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
- Drop reads below the 36 bases long (MINLEN:36)

### Step 2: Run `02_fastqc.sh`
Set input (first path) and output file location (-o) for forward (e.g. sample_R1_001_trimmed.fastq.gz) and reverse (e.g. sample_R2_001_trimmed.fastq.gz) reads.

### Step 3: Align genome with `03_salmon.sh`
Salmon is used for genome pseudo-alignment. Quant files are exported by Salmon and imported into R using tximport.

### Step 4: Perform Analysis in R using `04_analyze-bulk-seq-data.R` as a template.
An .Rmd script, `04_analyze-bulk-seq-data.Rmd` is also provided for differential expression at transcript- and gene-level analyses.

