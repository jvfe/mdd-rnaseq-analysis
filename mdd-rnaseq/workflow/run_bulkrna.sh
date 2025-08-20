#!/usr/bin/bash
set -e

# Run the pipeline with adjusted paths
nextflow run workflow/bulkrna \
    --input data/samplesheet.csv \
    --index data/references/index.idx \
    --tx2gene data/references/tx2gene.tsv \
    --gtf data/references/chr22.gtf \
    --outdir results \
    -profile docker \
    -resume
