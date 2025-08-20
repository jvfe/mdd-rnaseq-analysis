# README

## 1. Data sources

This task is based on publicly available sequencing data from the study "Sex-specific gene expression differences in the prefrontal cortex of major depressive disorder individuals". The study compares transcriptional changes in major depressive disorder and controls in both males and females. The dataset was originally sequenced using **\[Illumina HiSeq 2000]**.

The subsampled FASTQs are stored in `data/` and are used as the inputs for the workflow.

---

## 2. How to download

Data acquired from SRA

```bash
cd mdd-rnaseq/data
./download.sh
./download_references.sh
```

---

## 3. Pre-processing / subsampling

Filter down to only 12.5% of reads

```bash
cd mdd-rnaseq/data
./subsample.sh
```

---

## 4. How the workflow works

The workflow files is stored in `workflow/`.

---

### Step 1 – Re-run nextflow workflow described in article

**Purpose:** Remove low-quality reads and adapter sequences, align to a reference and generate a count matrix
**Tools:** `nextflow`, `fastqc`, `kallisto`, `tximport`, `multiqc`
**Inputs:** Subsampled FASTQ files (from`data/sra_data_downsampled/`), references (from `data/references/`)
**Outputs:** QC reports, Alignments, Count Matrix
**Command:**

```bash
bash workflow/run_bulkrna.sh
```

---

### Step 2 - Perform differential expression analysis and enrichment

**Purpose:** Find differentially expressed genes and transcripts between the conditions and perform functional enrichment
**Tools:** `R`, `edgeR`, `clusterProfiler`
**Inputs:** Count Matrix, metadata table
**Outputs:** Table of differentially expressed genes and trascripts
**Command:**

```bash
Rscript workflow/run_analysis.R
```
