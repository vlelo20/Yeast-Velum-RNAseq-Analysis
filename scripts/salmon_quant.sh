#!/usr/bin/env bash
set -euo pipefail

# Directories based on your 'tree' output
FASTQ_DIR="data/fastq"
INDEX="reference/salmon_index"
OUT_DIR="results/salmon_quant"
SAMPLE_SHEET="samplesheet.csv"

mkdir -p "$OUT_DIR"

echo "=== Starting Salmon Quantification ==="

# Read the samplesheet (skipping header)
tail -n +2 "$SAMPLE_SHEET" | while IFS=, read -r sample srr stage; do

    # Define the path to your single-end file
    # Based on your tree: data/fastq/IL20_1.fastq.gz
    READ_FILE="${FASTQ_DIR}/${sample}_1.fastq.gz"

    if [[ -f "$READ_FILE" ]]; then
        echo ">>> Processing Sample: ${sample} (${stage})"
        
        salmon quant -i "$INDEX" -l A \
            -r "$READ_FILE" \
            --validateMappings --gcBias --seqBias \
            -p 16 \
            -o "${OUT_DIR}/${sample}_${stage}"
    else
        echo "Error: File not found: $READ_FILE"
    fi

done

echo "=== All samples processed! ==="
