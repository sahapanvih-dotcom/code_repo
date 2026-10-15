#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# Advisor round 2, WP6: main Table 1 is generated from the candidate evidence master, so
# the shortlist and every statement about it regenerate from one release table.
#
# The locus, product, constraint percentile, triage class and the evidence tail (human
# homology, orthology grade, in-paralogues, C. albicans orthologue essentiality) are READ
# from output/table/rev_c2_candidate_evidence_master.csv. Only the curated rationale and the
# proposed first experiment are prose, and they are held here under version control so the
# table cannot drift from the data silently.
#
# Writes: output/table/rev_c3_main_table1.csv and build/_table1.md (markdown for the paper)
# Run: Rscript analysis_code/rev_c3_main_table1.R
# ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table)})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
m <- fread("output/table/rev_c2_candidate_evidence_master.csv")

ROWS <- list(
 AFUA_1G13300 = list(
  label = "GTP cyclohydrolase II, riboflavin",
  why = "Closest counterpart of *Candida* RIB1 (61% identity, reciprocal, conserved synteny across *Aspergillus*), but the gene tree places it and AFUA_2G01220 as **co-orthologues**, so the published *in vivo* validation attaches to the riboflavin pathway and not to this locus. Pathway absent from humans. Outside the ranked 815 on the prespecified power gate (14.7 expected variants, gate 15)",
  exp = "Double deletion with AFUA_2G01220; riboflavin auxotrophy"),
 AFUA_2G01220 = list(
  label = "GTP cyclohydrolase II, divergent paralogue",
  why = "Pipeline rank 12, Tier-1 (pLDDT 88, pocket 0.85); 29% identical to RIB1 over a third of its length, with an extra N-terminal domain and a *Pezizomycotina*-restricted two-gene cassette. Redundancy with the locus above is not excluded",
  exp = "As above — the pair, not either gene alone"),
 AFUA_8G05210 = list(
  label = "class I SAM-dependent methyltransferase",
  why = "Most constrained gene in the genome: 10 variants across 315 strains, none protein-changing. No functional hypothesis and no chemistry",
  exp = "Deletion with heterokaryon rescue"),
 AFUA_2G03560 = list(
  label = "two-component histidine kinase Bos1",
  why = "Constrained in both genera (*C. tropicalis* 83.6 within the orthologue-pair set); mammals have no two-component system, so selectivity is a class property",
  exp = "Deletion; survey agricultural pre-selection first"),
 AFUA_4G13510 = list(
  label = "isocitrate lyase 1",
  why = "Excluded in error by the paralogue filter, which matched it to methylisocitrate lyase at 46% identity — a different enzyme on a different substrate — and it is more constrained than malate synthase (Text S1 §S2.4). Outside the ranked 815 for that reason, and carried here with the correction stated",
  exp = "Glyoxylate flux under host-like carbon sources"),
 AFUA_4G10460 = list(
  label = "homocitrate synthase",
  why = "Both genera (*C. tropicalis* 95.5 within the pair set, 78% identity); α-aminoadipate lysine pathway, fungal-specific",
  exp = "¹³C/¹⁵N labelling before any deletion"),
 AFUA_6G03540 = list(
  label = "malate synthase",
  why = "Both genera (*C. tropicalis* 93.2 within the pair set); glyoxylate cycle absent from mammals; *M. tuberculosis* precedent",
  exp = "Glyoxylate-cycle flux under host-like conditions"),
 AFUA_1G02590 = list(
  label = "RNA polymerase I transcription factor",
  why = "Both genera (*C. tropicalis* 93.9 within the pair set); the strongest published external support in the intersection",
  exp = "Conditional depletion; assess human-counterpart risk first"),
 AFUA_1G09840 = list(
  label = "mitochondrial translocase TIM44",
  why = "Both genera (*C. tropicalis* 86.7 within the pair set)",
  exp = "Conditional depletion; assess human-counterpart risk first"),
 AFUA_5G09750 = list(
  label = "NCS2 family purine permease",
  why = "Tier-1, with homologues co-crystallised with purchasable purine analogues; no one-to-one *C. tropicalis* orthologue, and transporter redundancy is likely",
  exp = "Uptake assay"),
 AFUA_2G05060 = list(
  label = "alternative oxidase",
  why = "**Listed on genomic evidence only, outside the ranked 815.** Humans carry no alternative oxidase gene at all, so selectivity is a property of host biology rather than of compound design; no detectable human homologue, no *A. fumigatus* paralogue, present in all nine other *Aspergillus* species surveyed. It carries 12.6 expected variants, below the prespecified gate of 15, so it is not a ranked candidate and was not assessed for fold confidence or pocket druggability; it is not constrained in *C. tropicalis* (22.1 within the pair set, Table S10). Published work couples respiratory state to azole action in this organism [44] and a shelf inhibitor exists (SHAM [43]). **Drug-response induction of this gene was not measured in this study; the adjuvant hypothesis is prospective**",
  exp = "Prospective: qRT-PCR under voriconazole, then a voriconazole × SHAM checkerboard"),
 AFUA_3G10310 = list(
  label = "acetolactate synthase",
  why = "Listed so that it is not re-proposed: the strongest pocket in the set, but a crowded field with decades of agricultural pre-selection",
  exp = "**Do not start here** (Text S1 §S2.3)"))

TRIAGE_WORD <- c(eligible = "Eligible", conditional = "Conditional", exploratory = "Exploratory")

rows <- rbindlist(lapply(names(ROWS), function(L) {
  x <- m[locus_tag == L]
  stopifnot(nrow(x) == 1L)
  # evidence tail, built entirely from the master table
  tail <- character(0)
  if (x$homology_decision %like% "^cut")
    tail <- c(tail, sprintf("**cut by the high-sensitivity human screen** (%s, %.1f%% identity over %.0f%% of the protein, E = %.0e)",
                            sub("^sp\\|[^|]+\\|", "", x$homology_search_sensitive),
                            x$homology_sens_pident, x$homology_sens_qcov, x$homology_sens_evalue))
  tail <- c(tail, sprintf("orthology %s", x$orthology_confidence))
  if (x$n_af_inparalogs >= 1)
    tail <- c(tail, sprintf("%d in-paralogue(s) under permissive search", x$n_af_inparalogs))
  if (!is.na(x$n_src_reporting) && x$n_src_reporting > 0)
    tail <- c(tail, sprintf("*C. albicans* orthologue essential in %d of %d reporting sources",
                            x$n_src_essential, x$n_src_reporting))
  else
    tail <- c(tail, "no one-to-one *C. albicans* orthologue, so no external essentiality evidence")
  data.table(
    locus_tag = L,
    product_pctile = sprintf("%s (%.1f)", ROWS[[L]]$label, x$constraint_pctile),
    triage = TRIAGE_WORD[[x$triage]],
    evidence = { tl <- paste(tail, collapse = "; ")
                 # capitalise the first LETTER, which may sit inside markdown emphasis
                 i <- regexpr("[A-Za-z]", tl)
                 substr(tl, i, i) <- toupper(substr(tl, i, i))
                 paste0(ROWS[[L]]$why, ". ", tl) },
    first_experiment = ROWS[[L]]$exp,
    constraint_pctile = x$constraint_pctile, expected_variants = x$expected_variants,
    orthology_confidence = x$orthology_confidence, n_af_inparalogs = x$n_af_inparalogs,
    n_src_essential = x$n_src_essential, n_src_reporting = x$n_src_reporting,
    eligible_corrected = x$eligible_corrected, removal_reason = x$removal_reason)
}))
fwrite(rows, "output/table/rev_c3_main_table1.csv")

md <- c("| Locus | Product (Af constraint percentile) | Triage | Evidence and rationale | First experiment |",
        "|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s |", rows$locus_tag, rows$product_pctile,
                rows$triage, rows$evidence, rows$first_experiment))
writeLines(md, "doc/manuscript_2026-10/build/_table1.md")
cat(sprintf("wrote %d rows | eligible %d, conditional %d, exploratory %d\n", nrow(rows),
            sum(rows$triage == "Eligible"), sum(rows$triage == "Conditional"),
            sum(rows$triage == "Exploratory")))

## ---------------------------------------------------------------- splice into the paper
# build_asm.py lifts Table 1 out of the Results section of manuscript_main_v3.md, so the
# generated rows are written back into that file rather than kept beside it.
MS <- "doc/manuscript_2026-10/manuscript_main_v3.md"
src <- readLines(MS, warn = FALSE)
h <- grep("^\\| Locus \\| Product \\(Af constraint percentile\\)", src)
stopifnot(length(h) == 1L)
e <- h
while (e < length(src) && grepl("^\\|", src[e + 1L])) e <- e + 1L
writeLines(c(src[seq_len(h - 1L)], md, src[(e + 1L):length(src)]), MS)
cat(sprintf("spliced Table 1 into %s (replaced lines %d-%d)\n", MS, h, e))
