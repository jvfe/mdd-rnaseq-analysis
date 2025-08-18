#!/bin/bash

set -e

echo "--- Setting up directories ---"
mkdir -p sra_data reference_files

# -----------------------------------------------------------------------------
# 1. Download SRA Data using SRA Toolkit
# -----------------------------------------------------------------------------
# The file 'SRR_Acc_List.txt' should be in the same directory as this script.
# It contains the following accession numbers:

echo "--- Downloading SRA data ---"
cat SRR_Acc_List.txt | while read srr; do
  echo "Get $srr to FASTQ..."
  fastq-dump --split-files --gzip --outdir sra_data $srr
done

echo "--- SRA data download complete. Files are in sra_data/ ---"


