#!/usr/bin/env bash
# rev_b4_02c_domain_tree.sh
# Revision analysis B4, step 2c: the methodologically cleaner orthology test.
#
# WHY: the PF12471 (GTP_CH_N) N-terminal block carried by the AFUA_2G01220-type class
# has NO homologue in the Saccharomycotina or bacterial RIB1 proteins, so in the
# full-length alignment it is a large non-homologous block that inflates that class's
# branch (3.5-7.5 subs/site in the full-length trees - saturated, long-branch-attraction
# prone). Restricting the alignment to the PF00925 catalytic domain envelope - the only
# region homologous across ALL 29 family members - removes that confound.
# This is reported alongside, not instead of, the full-length result.
set -euo pipefail

BASE=/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta
OG=$BASE/analysis/orthogroup
LOG=$OG/logs
NT=${NT:-8}
PFAM=/mnt/HuangLab_26T/HuangLab/localDatabase/pfam/Pfam-A.hmm

# -- extract the PF00925 envelope (hmmsearch envelope coords, union of segments) --
hmmsearch --cpu $NT -E 1e-5 --domE 1e-3 \
  --domtblout "$OG"/pfam/rev_b4_family_PF00925.domtbl \
  "$OG"/pfam/PF00925.hmm "$OG"/fasta/rev_b4_family_labelled.faa \
  > "$LOG"/hmmsearch_family_PF00925.log 2>&1

awk '!/^#/ {print $1"\t"$20"\t"$21}' "$OG"/pfam/rev_b4_family_PF00925.domtbl \
  | sort -k1,1 -k2,2n \
  | awk -F'\t' '{if($1!=p){if(p!="")print p"\t"lo"\t"hi; p=$1; lo=$2; hi=$3} else {if($2<lo)lo=$2; if($3>hi)hi=$3}} END{if(p!="")print p"\t"lo"\t"hi}' \
  > "$OG"/pfam/rev_b4_PF00925_envelopes.tsv
cat "$OG"/pfam/rev_b4_PF00925_envelopes.tsv

awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2-1,$3,$1"_PF00925",0,"+"}' \
      "$OG"/pfam/rev_b4_PF00925_envelopes.tsv > "$OG"/pfam/rev_b4_PF00925_envelopes.bed
seqkit subseq --bed "$OG"/pfam/rev_b4_PF00925_envelopes.bed \
      "$OG"/fasta/rev_b4_family_labelled.faa 2> "$LOG"/seqkit_subseq.log \
  | sed -E 's/^>([^ ]+)_[0-9]+-[0-9]+:[+-].*$/>\1/' > "$OG"/fasta/rev_b4_family_PF00925dom.faa
seqkit fx2tab -nil "$OG"/fasta/rev_b4_family_PF00925dom.faa

mafft --localpair --maxiterate 1000 --thread $NT --anysymbol \
      "$OG"/fasta/rev_b4_family_PF00925dom.faa \
      > "$OG"/aln/rev_b4_PF00925dom_mafft_linsi.faa 2> "$LOG"/mafft_dom.log

trimal -in "$OG"/aln/rev_b4_PF00925dom_mafft_linsi.faa \
       -out "$OG"/trim/rev_b4_PF00925dom_automated1.faa -automated1 \
       > "$LOG"/trimal_dom.log 2>&1
printf "domain alignment ncol untrimmed=%s trimmed=%s\n" \
  "$(seqkit fx2tab -nl "$OG"/aln/rev_b4_PF00925dom_mafft_linsi.faa | head -1 | cut -f2)" \
  "$(seqkit fx2tab -nl "$OG"/trim/rev_b4_PF00925dom_automated1.faa | head -1 | cut -f2)"

iqtree2 -s "$OG"/trim/rev_b4_PF00925dom_automated1.faa -m MFP -B 1000 -bnni -alrt 1000 \
        -T $NT --seed 20261012 -o Ecol_NP4157931 \
        --prefix "$OG"/tree/rev_b4_domain -redo > "$OG"/tree/rev_b4_domain.stdout 2>&1
grep -m1 "Best-fit model" "$OG"/tree/rev_b4_domain.log
cat "$OG"/tree/rev_b4_domain.treefile
