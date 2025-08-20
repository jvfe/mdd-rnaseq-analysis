#!/bin/bash

# ==============================================================================
# Download and Subset Human Reference Genome from Ensembl
#
# This script performs the following actions:
# 1. Sets up variables for the latest Ensembl release URLs and output filenames.
# 2. Creates a directory to store the reference files.
# 3. Downloads the compressed FASTA (genome) and GTF (annotation) files.
# 4. Decompresses the downloaded files using a universally compatible method.
# 5. Subsets the GTF file to keep only annotations for chromosomes 1, 2, and 3.
# 6. Subsets the FASTA file using 'awk' to avoid dependency on 'samtools'.
# 7. Cleans up the large, intermediate files to save space.
#
# Usage:
#   - Save this script as 'download_reference.sh'
#   - Make it executable: chmod +x download_reference.sh
#   - Run it: ./download_reference.sh
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status.

FASTA_URL="http://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="http://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"

OUTPUT_DIR="references"
FASTA_FULL="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FULL="Homo_sapiens.GRCh38.112.gtf"
FASTA_SUBSET="Homo_sapiens.GRCh38.chr22.fa"
GTF_SUBSET="Homo_sapiens.GRCh38.112.chr22.gtf"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

if [ ! -f "${FASTA_FULL}.gz" ]; then
    wget -c "$FASTA_URL"
else
    echo "FASTA file already downloaded."
fi

if [ ! -f "${GTF_FULL}.gz" ]; then
    wget -c "$GTF_URL"
else
    echo "GTF file already downloaded."
fi

if [ ! -f "$FASTA_FULL" ]; then
    gunzip -c "${FASTA_FULL}.gz" > "$FASTA_FULL"
else
    echo "FASTA file already decompressed."
fi

if [ ! -f "$GTF_FULL" ]; then
    gunzip -c "${GTF_FULL}.gz" > "$GTF_FULL"
else
    echo "GTF file already decompressed."
fi

echo "--- Subsetting GTF file for chromosomes 22 ---"
# Keep header lines (starting with '#') and then filter the rest for chromosomes 1 in the first column.
(grep "^#" "$GTF_FULL"; awk -F'\t' 'BEGIN{OFS="\t"} $1 ~ /^[22]$/' "$GTF_FULL") > "$GTF_SUBSET"
echo "Subset GTF saved to: ${GTF_SUBSET}"

echo "--- Subsetting FASTA file for chromosomes 22 (using awk)... ---"
awk '
/^>/ {
    # Check if the header matches chromosome 22,
    # The match is for ">1 ", ">2 ", etc. to avoid matching ">10"
    if ($1 == ">22") {
        print_seq=1
    } else {
        print_seq=0
    }
}
# If the flag is set, print the line
print_seq {
    print
}' "$FASTA_FULL" > "$FASTA_SUBSET"

echo "Subset FASTA saved to: ${FASTA_SUBSET}"


echo "--- Cleaning up large intermediate files... ---"
rm "$FASTA_FULL"
rm "$GTF_FULL"
rm "${FASTA_FULL}.gz"
rm "${GTF_FULL}.gz"

echo "--- Script finished successfully! ---"
echo "Your subsetted files are in the '${OUTPUT_DIR}' directory."


