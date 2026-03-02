#!/usr/bin/env bash
set -euo pipefail

# ---- config ----
THREADS=16
ROOT=~/binf6110/assignment2
REF=$ROOT/reference
RAW=$ROOT/raw_fastq
OUT=$ROOT/results
SHEET=$ROOT/samplesheet.csv

mkdir -p "$RAW" "$REF" "$OUT"/{fastqc_raw,fastqc_trim,multiqc,salmon_quant} "$ROOT"/logs

# 1) References (skip if present)
cd "$REF"
[[ -s Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz ]] || \
  wget https://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
[[ -s Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz ]] || \
  wget https://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
[[ -s Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz ]] || \
  wget https://ftp.ensembl.org/pub/release-110/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz

# 2) Salmon decoy-aware index (skip if present)
if [[ ! -d salmon_index ]]; then
  zcat Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz \
    | grep '^>' | cut -d ' ' -f1 | sed 's/>//' > decoys.txt
  cat Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
      Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz \
    > gentrome.fa.gz
  salmon index -t gentrome.fa.gz -d decoys.txt -p "$THREADS" -i salmon_index
fi
cd "$ROOT"

# 3) ENA download of FASTQs
cut -d, -f2 "$SHEET" | tail -n +2 > runs.txt
> ena_links.tsv
while read acc; do
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&limit=0" \
  | tail -n +2 >> ena_links.tsv
done < runs.txt

while IFS=$'\t' read -r run urls md5s; do
  IFS=';' read -ra arr <<< "$urls"
  for u in "${arr[@]}"; do
    fn=$(basename "$u")
    [[ -s "$RAW/${run}_${fn}" ]] && continue
    wget -c "https://$u" -O "$RAW/${run}_${fn}"
  done
done < ena_links.tsv

# 4) QC
fastqc "$RAW"/*.fastq.gz -o "$OUT/fastqc_raw" -t "$THREADS"
multiqc "$OUT/fastqc_raw" -o "$OUT/multiqc"

# 5) Quantification
while IFS=, read -r sample srr stage; do
  [[ "$sample" == "sample" ]] && continue
  R1=$(ls "$RAW/${srr}_*1.fastq.gz" 2>/dev/null || true)
  R2=$(ls "$RAW/${srr}_*2.fastq.gz" 2>/dev/null || true)
  outdir="$OUT/salmon_quant/${sample}_${stage}"
  [[ -d "$outdir" ]] && continue
  if [[ -n "$R1" && -n "$R2" ]]; then
    salmon quant -i "$REF/salmon_index" -l A -1 "$R1" -2 "$R2" \
      --validateMappings --gcBias --seqBias \
      -p "$THREADS" -o "$outdir"
  else
    R=$(ls "$RAW/${srr}_*.fastq.gz" | head -n1)
    salmon quant -i "$REF/salmon_index" -l A -r "$R" \
      --validateMappings --gcBias --seqBias \
      -p "$THREADS" -o "$outdir"
  fi
done < "$SHEET"
