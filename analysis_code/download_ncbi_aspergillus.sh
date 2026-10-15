#!/usr/bin/env bash
# Bulk-download all NCBI assemblies for 11 Aspergillus species via NCBI datasets CLI.
# Usage: download_ncbi_aspergillus.sh <datasets_bin> <acc_list_dir> <raw_data_dir>
set -u
DS=$1; ACC=$2; RAW=$3
for sp in fumigatus lentulus flavus oryzae terreus niger tubingensis nidulans calidoustus candidus ochraceus; do
  out=$RAW/Aspergillus_$sp; mkdir -p "$out"; cd "$out" || exit 1
  n=$(wc -l < "$ACC/$sp.txt")
  echo "[$(date '+%F %T')] START $sp ($n accessions)"
  for try in 1 2 3; do
    "$DS" download genome accession --inputfile "$ACC/$sp.txt" \
      --include genome,gff3,protein,cds,rna,seq-report --dehydrated --no-progressbar \
      --filename "$sp.zip" && break
    echo "  dehydrated download failed (try $try)"; sleep 30
  done
  unzip -o -q "$sp.zip" -d dl && rm -f "$sp.zip"
  for try in 1 2 3; do
    "$DS" rehydrate --directory dl --max-workers 10 --no-progressbar && break
    echo "  rehydrate failed (try $try), retrying"; sleep 30
  done
  # flatten: one folder per accession
  for d in dl/ncbi_dataset/data/GC*; do [ -d "$d" ] && mv "$d" "$out/"; done
  mv dl/ncbi_dataset/data/*.jsonl dl/ncbi_dataset/data/*.json "$out/" 2>/dev/null
  mv dl/md5sum.txt "$out/ncbi_md5sum.txt" 2>/dev/null; rm -rf dl
  got=$(ls -d "$out"/GC* 2>/dev/null | wc -l)
  echo "[$(date '+%F %T')] DONE  $sp: $got/$n accession folders, $(du -sh "$out" | cut -f1)"
done
echo "[$(date '+%F %T')] ALL DONE"
