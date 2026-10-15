#!/usr/bin/env bash
# rev_b4_02_align_trim_tree.sh
# Revision analysis B4, step 2: align / trim / infer the ML gene tree for the
# PF00925 (GTP cyclohydrolase II, RIB1) family assembled in step 1.
#
# Documented choices
#   ALIGNER  mafft v7.505 (2022/Apr/10), L-INS-i: --localpair --maxiterate 1000
#            on FULL-LENGTH proteins (primary). A second alignment is built with
#            muscle 5.2 as an aligner-sensitivity control.
#   TRIMMING primary rule  : trimal -automated1  (heuristic selected for ML phylogeny)
#            sensitivity   : trimal -gt 0.8 -st 0.001  (strict, 80% occupancy)
#            sensitivity   : trimal -gt 0.3           (permissive)
#   TREE     iqtree2 2.0.7, -m MFP (ModelFinder Plus, automatic model selection),
#            -B 1000 (ultrafast bootstrap, with -bnni) and -alrt 1000 (SH-aLRT).
#   ROOTING  primary: Escherichia coli RibA (NP_415793.1) bacterial outgroup.
#            Sensitivity to rooting is assessed in step 3 (R/ape) with the ML tree
#            re-rooted at the outgroup, at the midpoint and at the Saccharomycotina clade.
set -euo pipefail

BASE=/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta
OG=$BASE/analysis/orthogroup
LOG=$OG/logs
NT=${NT:-8}

FAM=$OG/fasta/rev_b4_family.faa

# -------------------------------------------------- 0. tidy, informative FASTA labels
# IQ-TREE/newick-safe tip labels: <SpeciesTag>_<accession>
Rscript --vanilla - "$FAM" "$OG/proteomes/seq2species.tsv" \
        "$OG/fasta/rev_b4_family_labelled.faa" "$OG/fasta/rev_b4_tiplabels.tsv" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
fa  <- readLines(args[1])
map <- read.delim(args[2], header = FALSE, quote = "",
                  col.names = c("acc","tag","species","assembly","len"))
hdr <- grep("^>", fa)
acc <- sub("^>([^ ]+).*$", "\\1", fa[hdr])
m   <- map[match(acc, map$acc), ]
tip <- paste0(m$tag, "_", gsub("[^A-Za-z0-9]", "", acc))
fa[hdr] <- paste0(">", tip)
writeLines(fa, args[3])
write.table(data.frame(tip = tip, acc = acc, tag = m$tag, species = m$species,
                       assembly = m$assembly, aa_len = m$len),
            args[4], sep = "\t", quote = FALSE, row.names = FALSE)
cat("labelled", length(hdr), "sequences\n")
RS

# -------------------------------------------------- 1. alignment (versioned)
mafft --version 2> "$LOG"/mafft_version.txt || true
mafft --localpair --maxiterate 1000 --thread $NT --anysymbol \
      "$OG"/fasta/rev_b4_family_labelled.faa \
      > "$OG"/aln/rev_b4_family_mafft_linsi.faa 2> "$LOG"/mafft_linsi.log

muscle -version > "$LOG"/muscle_version.txt 2>&1 || true
muscle -align "$OG"/fasta/rev_b4_family_labelled.faa \
       -output "$OG"/aln/rev_b4_family_muscle.faa > "$LOG"/muscle.log 2>&1 || \
  echo "WARN muscle failed; aligner-sensitivity control unavailable" >&2

# -------------------------------------------------- 2. trimming (documented rules)
trimal --version > "$LOG"/trimal_version.txt 2>&1 || true
trimal -in "$OG"/aln/rev_b4_family_mafft_linsi.faa \
       -out "$OG"/trim/rev_b4_mafft_automated1.faa -automated1 \
       > "$LOG"/trimal_automated1.log 2>&1
trimal -in "$OG"/aln/rev_b4_family_mafft_linsi.faa \
       -out "$OG"/trim/rev_b4_mafft_gt08.faa -gt 0.8 -st 0.001 \
       > "$LOG"/trimal_gt08.log 2>&1
trimal -in "$OG"/aln/rev_b4_family_mafft_linsi.faa \
       -out "$OG"/trim/rev_b4_mafft_gt03.faa -gt 0.3 \
       > "$LOG"/trimal_gt03.log 2>&1
if [ -s "$OG"/aln/rev_b4_family_muscle.faa ]; then
  trimal -in "$OG"/aln/rev_b4_family_muscle.faa \
         -out "$OG"/trim/rev_b4_muscle_automated1.faa -automated1 \
         > "$LOG"/trimal_muscle_automated1.log 2>&1
fi

for f in "$OG"/aln/rev_b4_family_mafft_linsi.faa "$OG"/aln/rev_b4_family_muscle.faa \
         "$OG"/trim/rev_b4_mafft_automated1.faa "$OG"/trim/rev_b4_mafft_gt08.faa \
         "$OG"/trim/rev_b4_mafft_gt03.faa "$OG"/trim/rev_b4_muscle_automated1.faa; do
  [ -s "$f" ] && printf "%-55s ncol=%s\n" "$(basename "$f")" "$(seqkit fx2tab -nl "$f" | head -1 | cut -f2)"
done

# -------------------------------------------------- 3. per-sequence alignment coverage
Rscript --vanilla "$BASE"/analysis_code/rev_b4_02a_coverage.R

# -------------------------------------------------- 4. ML trees
run_iqtree () {  # $1 = trimmed alignment, $2 = prefix
  iqtree2 -s "$1" -m MFP -B 1000 -bnni -alrt 1000 -T $NT --seed 20261012 \
          -o Ecol_NP4157931 --prefix "$2" -redo > "${2}.stdout" 2>&1 || \
  iqtree2 -s "$1" -m MFP -B 1000 -bnni -alrt 1000 -T $NT --seed 20261012 \
          --prefix "$2" -redo > "${2}.stdout" 2>&1
}
run_iqtree "$OG"/trim/rev_b4_mafft_automated1.faa "$OG"/tree/rev_b4_main
run_iqtree "$OG"/trim/rev_b4_mafft_gt08.faa       "$OG"/tree/rev_b4_sens_gt08
run_iqtree "$OG"/trim/rev_b4_mafft_gt03.faa       "$OG"/tree/rev_b4_sens_gt03
[ -s "$OG"/trim/rev_b4_muscle_automated1.faa ] && \
  run_iqtree "$OG"/trim/rev_b4_muscle_automated1.faa "$OG"/tree/rev_b4_sens_muscle

grep -m1 "Best-fit model" "$OG"/tree/rev_b4_main.log || true

# -------------------------------------------------- 5. topology test of the two hypotheses
# H1 (manuscript): Candida/Sc RIB1 group with the AFUA_1G13300-type clade.
# H2 (alternative): Candida/Sc RIB1 group with the AFUA_2G01220-type clade.
# Constrained ML trees are inferred under each hypothesis with the same model, then
# compared against the unconstrained ML tree with the AU test (10000 RELL replicates).
Rscript --vanilla "$BASE"/analysis_code/rev_b4_02b_constraints.R

MODEL=$(grep -m1 "^Best-fit model" "$OG"/tree/rev_b4_main.log | sed 's/.*: //; s/ .*//')
for h in H1 H2; do
  iqtree2 -s "$OG"/trim/rev_b4_mafft_automated1.faa -m "$MODEL" \
          -g "$OG"/tree/rev_b4_constraint_${h}.tre -T $NT --seed 20261012 \
          --prefix "$OG"/tree/rev_b4_${h} -redo > "$OG"/tree/rev_b4_${h}.stdout 2>&1 || true
done
cat "$OG"/tree/rev_b4_main.treefile "$OG"/tree/rev_b4_H1.treefile \
    "$OG"/tree/rev_b4_H2.treefile 2>/dev/null > "$OG"/tree/rev_b4_hypotheses.trees
iqtree2 -s "$OG"/trim/rev_b4_mafft_automated1.faa -m "$MODEL" \
        -z "$OG"/tree/rev_b4_hypotheses.trees -n 0 -zb 10000 -au -T $NT \
        --seed 20261012 --prefix "$OG"/tree/rev_b4_autest -redo \
        > "$OG"/tree/rev_b4_autest.stdout 2>&1 || true

# -------------------------------------------------- versions
{
  echo "date	$(date '+%y:%m:%d, %H:%M:%S')"
  echo "mafft	$(mafft --version 2>&1 | head -1)"
  echo "muscle	$(muscle -version 2>&1 | head -1)"
  echo "trimal	$(trimal --version 2>&1 | sed -n '2p')"
  echo "iqtree2	$(iqtree2 --version 2>&1 | head -1)"
  echo "seqkit	$(seqkit version)"
  echo "R	$(Rscript -e 'cat(R.version.string)')"
  echo "ape	$(Rscript -e 'cat(as.character(packageVersion("ape")))')"
} > "$LOG"/rev_b4_step02_versions.txt
cat "$LOG"/rev_b4_step02_versions.txt
