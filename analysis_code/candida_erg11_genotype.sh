#!/usr/bin/env bash
# Call Erg11p protein variants per isolate by applying each VCF to the reference ERG11 CDS.
# C. tropicalis is DIPLOID, so bcftools consensus with -H I applies IUPAC codes for het sites;
# we call both a homozygous-only consensus and record het positions separately.
set -uo pipefail
C=/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta/analysis/candida
REG="NW_003020060.1:138041-139627"
BCF="docker run --rm -u $(id -u):$(id -g) -v /mnt/HuangLab_26T:/mnt/HuangLab_26T quay.io/biocontainers/bcftools:1.24--h118bc1c_2"
out=$C/erg11; mkdir -p $out
for v in $C/vcf/*.vcf.gz; do
  a=$(basename $v .vcf.gz)
  [ -s "$out/$a.fa" ] && continue
  $BCF bash -c "samtools faidx $C/ref/ref.fna $REG 2>/dev/null | bcftools consensus -s - $v 2>/dev/null" > $out/$a.fa 2>/dev/null
  [ -s "$out/$a.fa" ] || rm -f $out/$a.fa
done
echo "consensus sequences: $(ls $out/*.fa 2>/dev/null | wc -l)"
