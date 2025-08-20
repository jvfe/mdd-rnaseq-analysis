#!/bin/bash
set -e

echo "--- Setting up downsampled data directory ---"
mkdir -p sra_data_downsampled

SEED=1024
# Fraction of reads to keep
FRACTION=0.05

echo "--- Starting downsampling of FASTQ files ---"

for r1 in sra_data/*_1.fastq.gz; do
  r2="${r1/_1.fastq.gz/_2.fastq.gz}"
  
  base=$(basename "$r1" _1.fastq.gz)
  
  echo "Downsampling $base..."
  
  # Downsample R1 and R2 files
  seqtk sample -s $SEED "$r1" $FRACTION | gzip > sra_data_downsampled/${base}_1.fastq.gz
  seqtk sample -s $SEED "$r2" $FRACTION | gzip > sra_data_downsampled/${base}_2.fastq.gz
done

echo "--- Downsampling complete. Files are in sra_data_downsampled/ ---"
