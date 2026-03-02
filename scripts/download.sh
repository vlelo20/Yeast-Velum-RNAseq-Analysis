#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
SAMPLE_SHEET="samplesheet.csv"  # Ensure this file exists in the same directory
OUT_DIR="data/fastq"
THREADS=14                   # Adjust based on your CPU

mkdir -p "$OUT_DIR"

# Check if sample sheet exists
if [[ ! -f "$SAMPLE_SHEET" ]]; then
    echo "Error: $SAMPLE_SHEET not found."
    exit 1
fi

echo "=== Starting Download Process ==="

# Skip the header (tail -n +2) and loop through the CSV
# Using IFS=, to parse the comma-separated values
tail -n +2 "$SAMPLE_SHEET" | while IFS=, read -r sample srr stage; do

    # Check if files already exist to avoid redundant downloads
    if ls "${OUT_DIR}/${sample}"_*.fastq.gz >/dev/null 2>&1; then
        echo "--> Skipping ${sample} (${srr}); files already present."
        continue
    fi

    echo "=== Processing Sample: ${sample} (SRR: ${srr}, Stage: ${stage}) ==="

    # 1. Download
    # --include-technical is often needed for some SRA records, 
    # but --split-files is the standard for paired-end.
    fasterq-dump "${srr}" \
        --split-files \
        --outdir "$OUT_DIR" \
        --threads "$THREADS" \
        --temp .

    # 2. Rename SRR files to Sample Names and Gzip
    # This loop handles both single-end (one file) and paired-end (two files)
    for f in "${OUT_DIR}/${srr}"*.fastq; do
        if [[ -f "$f" ]]; then
            # Determine suffix (_1, _2, or empty)
            suffix=$(echo "$f" | sed "s/.*${srr}//")
            dest_file="${OUT_DIR}/${sample}${suffix}.gz"
            
            echo "  > Compressing and renaming to $(basename "$dest_file")"
            gzip -c "$f" > "$dest_file"
            rm "$f" # Remove the uncompressed original
        fi
    done

done

echo "=== All downloads and compression complete ==="
