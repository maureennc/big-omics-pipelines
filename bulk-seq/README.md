# Bulk-RNA-Sequencing
Bulk RNA Sequencing Pipeline for Harris lab naive vs. T. gondii-infected samples. Rivanna HPC used for bash and R scripting.

### Step 1: Run trimmomatic.slurm
Use Trimmomatic to trim out adapters from fastq files. Nextera adapters were used in example experiment. If unsure of which adapters should be trimmed, can run FastQC and read report.

Customizeable parameters used:
- Remove leading low quality or N bases (below quality 3) (LEADING:3)
- Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
- Drop reads below the 36 bases long (MINLEN:36)

### Step 2: Run fastqc.sh
Set input (first path) and output file location (-o) for forward (e.g. sample_R1_001_trimmed.fastq.gz) and reverse (e.g. sample_R2_001_trimmed.fastq.gz) reads.

### Step 3: Align genome with salmon.slurm
Salmon is used for genome pseudo-alignment. Quant files are exported by Salmon and imported into R using tximport.

### Step 4: Analysis in R
.Rmd script provided for differential expression at transcript- and gene-level analyses.
