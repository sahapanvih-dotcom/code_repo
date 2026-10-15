#!/usr/bin/env bash
# Per-isolate pipeline for C. tropicalis: depth-limited read download -> trim -> map -> variants + depth.
# Reads are taken as the first N records of each mate file via HTTP range requests, so pairing is
# preserved and we transfer ~25x coverage instead of the ~150x deposited.
set -uo pipefail
RUN=$1; ACC=$2
C=/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta/analysis/candida
REF=$C/ref/ref.fna
MM="$HOME/bin/micromamba run -n candida"
BCF="docker run --rm -u $(id -u):$(id -g) -v /mnt/HuangLab_26T:/mnt/HuangLab_26T -v /tmp/candida_work:/tmp/candida_work quay.io/biocontainers/bcftools:1.24--h118bc1c_2 bcftools"
export MAMBA_ROOT_PREFIX=$HOME/micromamba
NREC=9600000          # 2.4M read pairs -> ~45x theoretical, ~30x mapped on 14.6 Mb
BYTES=600000000       # range bytes per mate (generous; truncated gzip is fine)
T=/tmp/candida_work/$RUN; mkdir -p $T   # local NVMe, not NFS
[ -f "$C/vcf/$ACC.vcf.gz" ] && { echo "$ACC SKIP (done)"; exit 0; }

# --- resolve FTP paths
URLS=$(python3 -c "
import urllib.request,time,sys
r='$RUN'
u=f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={r}&result=read_run&fields=fastq_ftp&format=tsv'
out=''
for i in range(6):                      # ENA throttles metadata queries under load -> retry
    try:
        t=urllib.request.urlopen(u,timeout=90).read().decode().strip().splitlines()
        if len(t)>1 and t[1].count(chr(9))>=1:
            out=t[1].split(chr(9))[1]
            if out: break
    except Exception:
        pass
    time.sleep((i+1)*10)
print(out)
")
[ -z "$URLS" ] && { echo "$ACC FAIL no-ftp"; rm -rf $T; exit 1; }
U1=$(echo "$URLS" | cut -d';' -f1); U2=$(echo "$URLS" | cut -d';' -f2)
[ -z "$U2" ] || [ "$U1" = "$U2" ] && { echo "$ACC FAIL not-paired"; rm -rf $T; exit 1; }

# ENA refuses connections under heavy concurrency -> retry with backoff
fetch(){ local url="$1" out="$2" try
  for try in 1 2; do
    curl -sL --max-time 420 --retry 1 --retry-delay 5 -r 0-$BYTES "https://$url" 2>/dev/null | zcat 2>/dev/null | head -n $NREC > "$out"
    if [ "$(wc -l < "$out")" -ge 4000000 ]; then return 0; fi
    sleep $((try*10))
  done
  return 1
}
fetch "$U1" "$T/r1.fq"
fetch "$U2" "$T/r2.fq"

n1=$(( $(wc -l < $T/r1.fq) / 4 )); n2=$(( $(wc -l < $T/r2.fq) / 4 ))
n=$(( n1 < n2 ? n1 : n2 ))
[ "$n" -lt 1000000 ] && { echo "$ACC FAIL too-few-reads($n)"; rm -rf $T; exit 1; }
# truncate both to the same record count so mates stay in register
head -n $((n*4)) $T/r1.fq > $T/a1.fq; head -n $((n*4)) $T/r2.fq > $T/a2.fq; rm -f $T/r1.fq $T/r2.fq

$MM fastp -i $T/a1.fq -I $T/a2.fq -o $T/t1.fq -O $T/t2.fq -w 3 -q 20 -l 50 \
   --json $C/qc/$ACC.fastp.json --html /dev/null > /dev/null 2>&1
rm -f $T/a1.fq $T/a2.fq
$MM bwa mem -t 4 -R "@RG\tID:$ACC\tSM:$ACC\tPL:ILLUMINA" $REF $T/t1.fq $T/t2.fq 2>/dev/null \
  | $MM samtools sort -@ 2 -o $T/s.bam - 2>/dev/null
rm -f $T/t1.fq $T/t2.fq
$MM samtools index $T/s.bam
$MM mosdepth -n -t 2 --by $C/ref/genes.bed $C/depth/$ACC $T/s.bam 2>/dev/null
docker run --rm -u $(id -u):$(id -g) -v /mnt/HuangLab_26T:/mnt/HuangLab_26T -v /tmp/candida_work:/tmp/candida_work quay.io/biocontainers/bcftools:1.24--h118bc1c_2 bash -c "bcftools mpileup -f $REF -a AD,DP -q 20 -Q 20 -d 200 -Ou $T/s.bam 2>/dev/null | bcftools call -mv -Oz -o $C/vcf/$ACC.vcf.gz" 2>/dev/null
$BCF index $C/vcf/$ACC.vcf.gz 2>/dev/null
cov=$($MM samtools depth -a $T/s.bam 2>/dev/null | awk '{s+=$3} END{printf "%.1f", s/NR}')
nv=$($BCF view -H $C/vcf/$ACC.vcf.gz 2>/dev/null | wc -l)
if [ "${nv:-0}" -lt 1000 ]; then echo "$ACC FAIL empty-vcf(cov=${cov}x variants=$nv)"; rm -f $C/vcf/$ACC.vcf.gz*; rm -rf $T; exit 1; fi
echo "$ACC OK reads=$n cov=${cov}x variants=$nv"
rm -rf $T
