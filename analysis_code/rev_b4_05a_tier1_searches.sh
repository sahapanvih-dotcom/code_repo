#!/usr/bin/env bash
# =============================================================================
# rev_b4_05a_tier1_searches.sh
# Companion external-tool driver for rev_b4_05_tier1_orthogroup_check.R
# (revision analysis B4.05: orthogroup-aware orthology confidence, Tier-1 set)
#
# R is the primary analysis language for this project; this script exists only
# because the required steps are sequence-search steps for which the correct
# tools are blastp / hmmscan / seqkit (no adequate R implementation exists).
# ALL parsing, statistics, flagging and table writing happen in the R script.
#
# Usage: bash rev_b4_05a_tier1_searches.sh <PROJECT_BASE>
# Expects <BASE>/analysis/orthogroup/tier1/fasta/tier1_af_accessions.txt
#   (one A. fumigatus RefSeq protein accession per line, written by the R script)
#
# Writes: <BASE>/analysis/orthogroup/tier1/{fasta,db,blast,pfam,logs}/
# =============================================================================
set -euo pipefail

BASE="${1:-/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta}"
THREADS="${THREADS:-32}"

AF_FAA="$BASE/analysis/target_discovery/hom/af293.faa"
CT_FAA="$BASE/analysis/candida/ref/protein.faa"
CALB_FAA="$BASE/analysis/orthogroup/ref/calb_sc5314.faa"
PFAM="/mnt/HuangLab_26T/HuangLab/localDatabase/pfam/Pfam-A.hmm"

WD="$BASE/${REV_B4_WORKDIR:-analysis/orthogroup/tier1}"   # overridable so the same audit can run over a wider candidate set
FA="$WD/fasta"; DB="$WD/db"; BL="$WD/blast"; PF="$WD/pfam"; LG="$WD/logs"
mkdir -p "$FA" "$DB" "$BL" "$PF" "$LG"

ACC="$FA/tier1_af_accessions.txt"
sync 2>/dev/null || true
for _t in 1 2 3 4 5 6 7 8 9 10; do [[ -s "$ACC" ]] && break; sleep 1; done
[[ -s "$ACC" ]] || { echo "ERROR: missing $ACC" >&2; exit 1; }

for f in "$AF_FAA" "$CT_FAA" "$CALB_FAA" "$PFAM"; do
  [[ -s "$f" ]] || { echo "ERROR: missing input $f" >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# IMPORTANT: blastp uses -evalue BOTH as the reporting threshold AND to derive the
# preliminary (ungapped) score cutoff. Running with -evalue 1e-5 therefore DROPS
# alignments whose final E-value is far below 1e-5 (verified here: XP_749067.1 vs
# C. tropicalis XP_002548511.1 scores E=1.45e-09 but is reported only when
# -evalue >= 1e-3). We therefore report permissively and apply the E<=1e-5
# criterion downstream in R.
# ---------------------------------------------------------------------------
REPORT_EVALUE="${REPORT_EVALUE:-10}"
EVAL_MAX="1e-5"   # significance cut applied in R (and for best-hit selection below)
OFMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs'

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LG/rev_b4_05a_run.log"; }

log "START rev_b4_05a_tier1_searches.sh  BASE=$BASE THREADS=$THREADS"

# ---------------------------------------------------------------- 1. candidates
log "STEP 1 seqkit grep: extract Tier-1 candidate proteins from af293.faa"
# NOTE: --quiet is essential - seqkit writes its "N patterns loaded" INFO line to
# STDOUT, which would otherwise be prepended to the FASTA and break hmmscan.
seqkit grep --quiet -f "$ACC" "$AF_FAA" > "$FA/tier1_candidates.faa" 2> "$LG/seqkit_candidates.log"
sync 2>/dev/null || true
n_got=$(grep -c '^>' "$FA/tier1_candidates.faa" || true)
n_want=$(grep -c . "$ACC")
log "  extracted $n_got / $n_want requested"
[[ "$n_got" -eq "$n_want" ]] || { echo "ERROR: seqkit extracted $n_got of $n_want candidates" >&2; exit 1; }

# ------------------------------------------------------------------ 2. BLAST DBs
log "STEP 2 makeblastdb (af293, C. tropicalis, C. albicans SC5314)"
for pair in "af293:$AF_FAA" "ctrop:$CT_FAA" "calb:$CALB_FAA"; do
  name="${pair%%:*}"; faa="${pair#*:}"
  if [[ ! -s "$DB/$name.pin" && ! -s "$DB/$name.00.pin" ]]; then
    makeblastdb -in "$faa" -dbtype prot -out "$DB/$name" -title "$name" \
      >> "$LG/makeblastdb.log" 2>&1
  fi
done

# ------------------------------------------- 3. forward: candidate -> af293 self
# (in-paralog detection; keep deep hit list so paralogue counting is not capped)
log "STEP 3 blastp candidate -> af293 (in-paralogs)"
blastp -query "$FA/tier1_candidates.faa" -db "$DB/af293" \
  -evalue "$REPORT_EVALUE" -outfmt "$OFMT" -max_target_seqs 2000 -num_threads "$THREADS" \
  -out "$BL/cand_vs_af293.tsv" 2>> "$LG/blastp.log"

# -------------------------------------------- 4. forward: candidate -> Candida x2
log "STEP 4 blastp candidate -> C. tropicalis / C. albicans"
blastp -query "$FA/tier1_candidates.faa" -db "$DB/ctrop" \
  -evalue "$REPORT_EVALUE" -outfmt "$OFMT" -max_target_seqs 500 -num_threads "$THREADS" \
  -out "$BL/cand_vs_ctrop.tsv" 2>> "$LG/blastp.log"
blastp -query "$FA/tier1_candidates.faa" -db "$DB/calb" \
  -evalue "$REPORT_EVALUE" -outfmt "$OFMT" -max_target_seqs 500 -num_threads "$THREADS" \
  -out "$BL/cand_vs_calb.tsv" 2>> "$LG/blastp.log"

# --------------------------------------------------- 5. reciprocal (reverse) BLAST
# top-1 subject per query by bitscore, then blastp that protein back to af293
log "STEP 5 reciprocal blastp (Candida best hit -> af293)"
for sp in ctrop calb; do
  src_faa=$CT_FAA; [[ "$sp" == "calb" ]] && src_faa=$CALB_FAA
  awk -F'\t' -v e="$EVAL_MAX" '$11+0 <= e+0' "$BL/cand_vs_${sp}.tsv" \
    | sort -k1,1 -k12,12gr | awk -F'\t' '!seen[$1]++ {print $2}' \
    | sort -u > "$FA/${sp}_besthits.ids"
  sync 2>/dev/null || true
  if [[ -s "$FA/${sp}_besthits.ids" ]]; then
    seqkit grep --quiet -f "$FA/${sp}_besthits.ids" "$src_faa" > "$FA/${sp}_besthits.faa" \
      2>> "$LG/seqkit_besthits.log"
    sync 2>/dev/null || true
    b_got=$(grep -c '^>' "$FA/${sp}_besthits.faa" || true)
    b_want=$(grep -c . "$FA/${sp}_besthits.ids")
    [[ "$b_got" -eq "$b_want" ]] || { echo "ERROR: $sp seqkit got $b_got of $b_want best hits" >&2; exit 1; }
    blastp -query "$FA/${sp}_besthits.faa" -db "$DB/af293" \
      -evalue "$REPORT_EVALUE" -outfmt "$OFMT" -max_target_seqs 500 -num_threads "$THREADS" \
      -out "$BL/${sp}_besthits_vs_af293.tsv" 2>> "$LG/blastp.log"
  else
    : > "$FA/${sp}_besthits.faa"; : > "$BL/${sp}_besthits_vs_af293.tsv"
  fi
  log "  $sp: $(wc -l < "$FA/${sp}_besthits.ids") distinct best hits"
done

# ------------------------------------------------------------ 6. Pfam architecture
sync 2>/dev/null || true
log "STEP 6 hmmscan vs Pfam-A (candidates and C. tropicalis best hits)"
run_hmmscan () {  # $1=query faa  $2=output tag
  if [[ -s "$PF/$2.domtbl" ]]; then log "  $2.domtbl exists, skipping"; return; fi
  # reported permissively (HMMER default E=10); the domain i-Evalue<=1e-5 cut is
  # applied in R, for symmetry with the BLAST handling above
  hmmscan --cpu "$THREADS" \
    --domtblout "$PF/$2.domtbl" --tblout "$PF/$2.tbl" \
    "$PFAM" "$1" > "$PF/$2.hmmscan.out" 2>> "$LG/hmmscan.log"
}
run_hmmscan "$FA/tier1_candidates.faa" "tier1_candidates_pfam"
[[ -s "$FA/ctrop_besthits.faa" ]] && run_hmmscan "$FA/ctrop_besthits.faa" "ctrop_besthits_pfam"

# --------------------------------------------------------------- 7. tool versions
log "STEP 7 record tool versions"
{
  echo "# rev_b4_05 Tier-1 orthogroup check - exact tool version strings"
  echo "# generated $(date '+%Y-%m-%d %H:%M:%S %Z') on $(hostname)"
  echo
  echo "## blastp"            ; blastp -version 2>&1
  echo "## makeblastdb"       ; makeblastdb -version 2>&1
  echo "## hmmscan"           ; hmmscan -h 2>&1 | sed -n '2p'
  echo "## seqkit"            ; seqkit version 2>&1
  echo "## diamond (available, not used)"; diamond --version 2>&1
  echo "## R"                 ; Rscript --version 2>&1
  echo "## R packages"        ; Rscript -e 'for (p in c("ape")) cat(p, as.character(packageVersion(p)), "\n")' 2>/dev/null || true
  echo
  echo "## databases"
  echo "Pfam-A.hmm           : $PFAM  ($(stat -c%s "$PFAM") bytes, mtime $(stat -c%y "$PFAM"))"
  echo "af293.faa            : $AF_FAA ($(grep -c '^>' "$AF_FAA") seqs)"
  echo "C. tropicalis faa    : $CT_FAA ($(grep -c '^>' "$CT_FAA") seqs)"
  echo "C. albicans SC5314   : $CALB_FAA ($(grep -c '^>' "$CALB_FAA") seqs)"
  echo
  echo "## exact command lines run"
  echo "seqkit grep --quiet -f tier1_af_accessions.txt af293.faa > tier1_candidates.faa"
  echo "makeblastdb -in <faa> -dbtype prot -out <db> -title <db>"
  echo "blastp -query tier1_candidates.faa -db af293 -evalue $REPORT_EVALUE -outfmt '$OFMT' -max_target_seqs 2000 -num_threads $THREADS -out cand_vs_af293.tsv"
  echo "blastp -query tier1_candidates.faa -db ctrop -evalue $REPORT_EVALUE -outfmt '$OFMT' -max_target_seqs 500 -num_threads $THREADS -out cand_vs_ctrop.tsv"
  echo "blastp -query tier1_candidates.faa -db calb  -evalue $REPORT_EVALUE -outfmt '$OFMT' -max_target_seqs 500 -num_threads $THREADS -out cand_vs_calb.tsv"
  echo "awk -F'\\t' '\$11+0 <= $EVAL_MAX' cand_vs_<sp>.tsv | sort -k1,1 -k12,12gr | awk '!seen[\$1]++ {print \$2}' | sort -u > <sp>_besthits.ids"
  echo "seqkit grep --quiet -f <sp>_besthits.ids <sp>.faa > <sp>_besthits.faa"
  echo "blastp -query <sp>_besthits.faa  -db af293 -evalue $REPORT_EVALUE -outfmt '$OFMT' -max_target_seqs 500 -num_threads $THREADS -out <sp>_besthits_vs_af293.tsv"
  echo "hmmscan --cpu $THREADS --domtblout <tag>.domtbl --tblout <tag>.tbl Pfam-A.hmm <query.faa>"
  echo
  echo "## significance criteria (applied in rev_b4_05_tier1_orthogroup_check.R, not by the tools)"
  echo "BLAST: E <= $EVAL_MAX ; in-paralog additionally requires qcovs >= 50"
  echo "Pfam : per-domain i-Evalue <= $EVAL_MAX, overlapping domains resolved best-score-first"
  echo "Rationale: blastp -evalue sets the preliminary extension cutoff as well as the"
  echo "reporting threshold, so searching at -evalue 1e-5 silently loses alignments whose"
  echo "final E-value is below 1e-5 (e.g. XP_749067.1 vs XP_002548511.1, E=1.45e-09)."
} > "$WD/rev_b4_tier1_versions.txt"

sync 2>/dev/null || true
log "DONE rev_b4_05a_tier1_searches.sh"
