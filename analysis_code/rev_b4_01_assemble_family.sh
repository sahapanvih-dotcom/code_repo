#!/usr/bin/env bash
# rev_b4_01_assemble_family.sh
# Revision analysis B4, step 1: assemble the complete GTP cyclohydrolase II (RIB1)
# protein family across 10 annotated Aspergillus species, C. tropicalis, C. albicans,
# S. cerevisiae (fungal outgroup) and E. coli (bacterial outgroup / root).
#
# External tools are used where they are the right tool for the step (blastp, hmmscan,
# seqkit); R is used downstream for statistics, summarising and plotting.
set -euo pipefail

BASE=/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta
OG=$BASE/analysis/orthogroup
RAW=$BASE/raw_data/coverage_pass
LOG=$OG/logs
mkdir -p "$OG"/{fasta,aln,trim,tree,logs,pfam,synteny,blast,ref,proteomes}

PFAM=/mnt/HuangLab_26T/HuangLab/localDatabase/pfam/Pfam-A.hmm

# ---------------------------------------------------------------- 1. proteome set
# One annotated representative assembly per Aspergillus species (RefSeq GCF preferred).
cat > "$OG"/proteomes/representatives.tsv <<'EOF'
Aspergillus_calidoustus	GCA_001511075.1	Acal
Aspergillus_candidus	GCF_002847045.1	Acan
Aspergillus_flavus	GCF_009017415.1	Afla
Aspergillus_fumigatus	GCF_000002655.1	Afum
Aspergillus_lentulus	GCF_010724455.1	Alen
Aspergillus_nidulans	GCF_000011425.1	Anid
Aspergillus_niger	GCF_000002855.4	Anig
Aspergillus_oryzae	GCF_000184455.2	Aory
Aspergillus_terreus	GCF_000149615.1	Ater
Aspergillus_tubingensis	GCF_013340325.1	Atub
EOF
# NOTE: Aspergillus_ochraceus is in raw_data/ but none of its 4 assemblies carries
# protein.faa / genomic.gff, so it cannot contribute an annotated protein; it is handled
# separately by tblastn in step 1c (evidence only, not placed in the tree).

: > "$OG"/proteomes/searchdb.faa
: > "$OG"/proteomes/seq2species.tsv
while IFS=$'\t' read -r sp acc tag; do
  faa="$RAW/$sp/$acc/protein.faa"
  [ -f "$faa" ] || { echo "MISSING $faa" >&2; exit 1; }
  seqkit fx2tab -nil "$faa" | awk -v t="$tag" -v s="$sp" -v a="$acc" -F'\t' \
    '{split($1,h," "); print h[1]"\t"t"\t"s"\t"a"\t"$2}' >> "$OG"/proteomes/seq2species.tsv
  cat "$faa" >> "$OG"/proteomes/searchdb.faa
done < "$OG"/proteomes/representatives.tsv

# Candida tropicalis MYA-3404 (the reference behind the CTRG_ ortholog map)
for spec in "Ctro:Candida_tropicalis:GCF_000006335.3:$BASE/analysis/candida/ref/protein.faa" \
            "Calb:Candida_albicans:GCF_000182965.3:$OG/ref/calb_sc5314.faa" \
            "Scer:Saccharomyces_cerevisiae:GCF_000146045.2:$OG/ref/scer_s288c.faa" \
            "Ecol:Escherichia_coli:GCF_000005845.2:$OG/ref/ecoli_k12.faa"; do
  IFS=: read -r tag sp acc faa <<<"$spec"
  [ -f "$faa" ] || { echo "MISSING $faa" >&2; exit 1; }
  seqkit fx2tab -nil "$faa" | awk -v t="$tag" -v s="$sp" -v a="$acc" -F'\t' \
    '{split($1,h," "); print h[1]"\t"t"\t"s"\t"a"\t"$2}' >> "$OG"/proteomes/seq2species.tsv
  cat "$faa" >> "$OG"/proteomes/searchdb.faa
done

echo "searchdb sequences: $(grep -c '>' "$OG"/proteomes/searchdb.faa)"

# ---------------------------------------------------------------- 2. seeds
# Both A. fumigatus candidates + the published Candida RIB1 proteins + Sc/Ec anchors.
cat > "$OG"/fasta/rev_b4_seed_ids.txt <<'EOF'
XP_077660111.1
XP_749258.2
XP_002546018.1
XP_715792.1
NP_009520.1
NP_415793.1
EOF
seqkit grep -r -f "$OG"/fasta/rev_b4_seed_ids.txt "$OG"/proteomes/searchdb.faa \
  > "$OG"/fasta/rev_b4_seeds.faa
echo "seeds: $(grep -c '>' "$OG"/fasta/rev_b4_seeds.faa)"

# ---------------------------------------------------------------- 3. homology search
makeblastdb -in "$OG"/proteomes/searchdb.faa -dbtype prot \
  -out "$OG"/blast/searchdb -logfile "$LOG"/makeblastdb.log

blastp -query "$OG"/fasta/rev_b4_seeds.faa -db "$OG"/blast/searchdb \
  -evalue 1e-5 -max_target_seqs 2000 -num_threads 8 \
  -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
  > "$OG"/blast/rev_b4_seeds_vs_searchdb.tsv 2> "$LOG"/blastp_seeds.log
echo "seed blastp hits: $(wc -l < "$OG"/blast/rev_b4_seeds_vs_searchdb.tsv)"

# Candidate family members: any hit to any seed with E<=1e-5 covering >=40% of the seed.
awk -F'\t' '{qcov=($8-$7+1)/$5} $11<=1e-5 && qcov>=0.40 {print $2}' \
  "$OG"/blast/rev_b4_seeds_vs_searchdb.tsv | sort -u > "$OG"/fasta/rev_b4_candidates.ids
echo "blast candidates: $(wc -l < "$OG"/fasta/rev_b4_candidates.ids)"

seqkit grep -f "$OG"/fasta/rev_b4_candidates.ids "$OG"/proteomes/searchdb.faa \
  > "$OG"/fasta/rev_b4_candidates.faa

# ---------------------------------------------------------------- 4. Pfam domain check
# Pfam-A is available locally, so domain architecture is determined by hmmscan rather
# than inferred from annotation text.
hmmscan --cpu 8 -E 1e-5 --domE 1e-5 \
  --tblout "$OG"/pfam/rev_b4_candidates_pfam.tbl \
  --domtblout "$OG"/pfam/rev_b4_candidates_pfam.domtbl \
  "$PFAM" "$OG"/fasta/rev_b4_candidates.faa > "$LOG"/hmmscan_candidates.log 2>&1
echo "hmmscan done"

# Family definition: must carry the GTP_cyclohydro2 (PF00925) domain.
awk '!/^#/ && $1=="GTP_cyclohydro2" {print $3}' "$OG"/pfam/rev_b4_candidates_pfam.tbl \
  | sort -u > "$OG"/fasta/rev_b4_family.ids
echo "PF00925-confirmed family members: $(wc -l < "$OG"/fasta/rev_b4_family.ids)"

seqkit grep -f "$OG"/fasta/rev_b4_family.ids "$OG"/proteomes/searchdb.faa \
  > "$OG"/fasta/rev_b4_family.faa
seqkit fx2tab -nil "$OG"/fasta/rev_b4_family.faa > "$OG"/fasta/rev_b4_family_lengths.tsv

# ---------------------------------------------------------------- 1c. A. ochraceus tblastn
# A. ochraceus has genome but no protein annotation: record tblastn evidence for how many
# GTP cyclohydrolase II loci it carries (used as independent copy-number evidence only).
AO=$(ls -d "$RAW"/Aspergillus_ochraceus/*/ 2>/dev/null | head -1)
if [ -n "${AO:-}" ]; then
  G=$(ls "$AO"/*_genomic.fna 2>/dev/null | head -1)
  if [ -n "${G:-}" ]; then
    tblastn -query "$OG"/fasta/rev_b4_seeds.faa -subject "$G" -evalue 1e-10 \
      -outfmt '6 qseqid sseqid pident length qlen evalue bitscore sstart send' \
      > "$OG"/blast/rev_b4_ochraceus_tblastn.tsv 2> "$LOG"/tblastn_ochraceus.log || true
    echo "A. ochraceus tblastn hits: $(wc -l < "$OG"/blast/rev_b4_ochraceus_tblastn.tsv)"
  fi
fi

# ---------------------------------------------------------------- 5. host-homology check
blastp -query "$OG"/fasta/rev_b4_seeds.faa \
  -subject "$BASE/analysis/target_discovery/hom/human.faa" \
  -evalue 1e-3 -num_threads 4 \
  -outfmt '6 qseqid sseqid pident length qlen evalue bitscore' \
  > "$OG"/blast/rev_b4_seeds_vs_human.tsv 2> "$LOG"/blastp_human.log || true
echo "seed-vs-human hits (E<=1e-3): $(wc -l < "$OG"/blast/rev_b4_seeds_vs_human.tsv)"

# ---------------------------------------------------------------- versions
{
  echo "date	$(date '+%y:%m:%d, %H:%M:%S')"
  echo "blastp	$(blastp -version | head -1)"
  echo "makeblastdb	$(makeblastdb -version | head -1)"
  echo "hmmscan	$(hmmscan -h | sed -n '2p')"
  echo "seqkit	$(seqkit version)"
  echo "pfam_db	$PFAM"
} > "$LOG"/rev_b4_step01_versions.txt
cat "$LOG"/rev_b4_step01_versions.txt
