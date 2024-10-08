#!/bin/bash

# Define variables
FASTQ_PATH={'azenta/fastq_path'}

# Run cellranger count
cellranger count \
  --id=1 \
  --sample=1 \
  --fastqs={fastq_path} \
  --transcriptome=refdata-gex-mm10-2020-A \
  --localmem=80 \
  --localcores=16 \
  --r1-length=28 \
  --r2-length=98 \
  --chemistry=threeprime \
  --project=30-996613217
  
## Give the script executable permissions by running: `chmod +x run_cellranger.sh`
## All samples were sequenced on a single Nova X 25B lane