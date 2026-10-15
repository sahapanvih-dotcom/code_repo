#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# Advisor round 2, WP6 (with WP2): one release-ready candidate-level evidence table.
#
# Every statement the manuscript makes about a named candidate must be regenerable from
# this one table. It covers the union of the 18 Tier-1 candidates and the 12 genes carried
# in main Table 1 (26 genes), and reports the seven evidence domains the brief specifies:
# constraint, robustness, human homology, paralogy/orthology, structure, external
# validation, final triage.
#
# WP2 asks for non-overlapping triage labels under prespecified rules. The cascade below is
# evaluated in order and the first match wins, so the classes cannot overlap:
#
#   1 exploratory  - a human homologue is detectable at high sensitivity. This is the
#                    pipeline's own prespecified exclusion criterion, not a safety finding;
#                    such a gene is reported as a documented exclusion, never as a lead.
#   2 conditional  - fails a non-selectivity criterion for a reason stated in the
#                    manuscript (a documented exception, listed in DOCUMENTED_EXCEPTIONS).
#   3 exploratory  - fails any other criterion of the corrected funnel.
#   4 conditional  - passes the corrected funnel but carries orthology or paralogy risk:
#                    orthology confidence LOW or UNRESOLVED, or an in-paralogue detected
#                    under permissive search.
#   5 eligible     - passes the corrected funnel with orthology confidence HIGH or MEDIUM
#                    and no in-paralogue.
#
# Run: Rscript analysis_code/rev_c2_candidate_evidence_master.R
#      (needs rev_c1_canonical_filter_audit.R and the _union26 orthology audit first)
# ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table)})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
TAB <- "output/table"
say <- function(...) cat(sprintf(...), "\n", sep = "")

# Genes carried in main Table 1 that the funnel drops for a reason the manuscript states.
DOCUMENTED_EXCEPTIONS <- list(
  AFUA_1G13300 = "below the >=15 expected-variant power gate; carried because it is the counterpart of the published Candida RIB1 validation",
  AFUA_2G05060 = "below the >=15 expected-variant power gate; carried because humans have no alternative oxidase, so selectivity is structural rather than statistical",
  AFUA_4G13510 = "removed by the paralogue filter, which the cross-genus analysis showed to be wrong here; carried with that correction stated")

aud  <- fread(file.path(TAB, "rev_c1_canonical_filter_audit.csv"))
uni  <- fread(file.path(TAB, "rev_c2_candidate_set_union.csv"))
orth <- fread(file.path(TAB, "rev_b4_tier1_orthology_confidence_union26.csv"))
para <- fread(file.path(TAB, "rev_b4_tier1_paralog_detail_union26.csv"))
pct  <- fread(file.path(TAB, "rev_b2_candidate_percentiles.csv"))
ret  <- fread(file.path(TAB, "rev_b2_tier1_retention.csv"))
str_ <- fread("doc/Aspergillus/target_discovery_2026-08/tables/candidates_full_assessment.csv")
xg   <- fread(file.path(TAB, "ms_table7_crossgenus_candidates.csv"))

## ---------------------------------------------------------------- external validation
# From rev_c0_calb_essentiality_per_locus.R, which reproduces the identifier crosswalk and
# per-source admissibility rules of rev_b1_essentiality_benchmark.R. The published
# ms_table10_intersection_essentiality.csv predates the orf19 crosswalk and under-counts
# GRACEv1_OMeara2015; rev_c0 writes the agreement check.
esum <- fread(file.path(TAB, "rev_c0_calb_essentiality_per_af_locus.csv"))
setnames(esum, "af_locus", "locus_tag")
esum <- esum[, .(locus_tag, calb_acc, calb_pident, calb_qcov,
                 n_src_reporting = n_src_reporting_strict,
                 n_src_essential = n_src_essential_strict,
                 sources_essential)]

## ---------------------------------------------------------------- assemble
m <- merge(uni[, .(locus_tag, set, rank_published = fifelse(rank >= 9000L, NA_integer_, rank))],
           aud, by = "locus_tag", all.x = TRUE)
m <- merge(m, orth[, .(locus_tag, orthology_confidence, n_risk_flags, risk_flags,
                       n_af_inparalogs, best_inparalog_acc, best_inparalog_pident,
                       ct_locus_recomputed, ct_reciprocal, calb_locus_recomputed,
                       calb_reciprocal, domain_mismatch_class, confidence_rule)],
           by = "locus_tag", all.x = TRUE)
m <- merge(m, pct[, .(locus_tag, pctile_median, pctile_q05, pctile_q95,
                      rank_median, rank_q05, rank_q95, eligible_freq)],
           by = "locus_tag", all.x = TRUE)
m <- merge(m, str_[, .(locus_tag, mean_plddt, n_pockets, best_druggability,
                       best_pocket_volume, best_pdb_id, best_pdb, top_ligand)],
           by = "locus_tag", all.x = TRUE)
m <- merge(m, esum, by = "locus_tag", all.x = TRUE)
m <- merge(m, xg[, .(locus_tag = Af_locus, ct_orthologue = Ct_orthologue,
                     ct_identity_pct = identity_pct, ct_piNpiS = Ct_piNpiS,
                     ct_pctile = Ct_pctile, crossgenus_quadrant = Quadrant)],
           by = "locus_tag", all.x = TRUE)

## ---------------------------------------------------------------- triage cascade
m[, documented_exception := locus_tag %in% names(DOCUMENTED_EXCEPTIONS)]
m[, exception_reason := vapply(locus_tag, function(L)
    if (!is.null(DOCUMENTED_EXCEPTIONS[[L]])) DOCUMENTED_EXCEPTIONS[[L]] else "", character(1))]
m[, orthology_risk := orthology_confidence %in% c("LOW", "UNRESOLVED") | n_af_inparalogs >= 1]

m[, triage := fcase(
    !c4_no_human_sen,                        "exploratory",
    !eligible_corrected & documented_exception, "conditional",
    !eligible_corrected,                     "exploratory",
    orthology_risk,                          "conditional",
    default =                                "eligible")]
m[, triage_rationale := fcase(
    !c4_no_human_sen, sprintf("human homologue detectable at high sensitivity (%s, %.1f%% id, %.0f%% qcov, E=%.1e); prespecified exclusion criterion of this pipeline, not a toxicity measurement",
                              hs_sen_hit, hs_sen_pident, hs_sen_qcov, hs_sen_evalue),
    !eligible_corrected & documented_exception, exception_reason,
    !eligible_corrected, sprintf("fails the corrected funnel at %s with no stated exception", removal_reason),
    orthology_risk, sprintf("passes the corrected funnel; orthology confidence %s%s",
                            orthology_confidence,
                            fifelse(n_af_inparalogs >= 1,
                                    sprintf(", %d in-paralogue(s) detected under permissive search", n_af_inparalogs), "")),
    default = sprintf("passes the corrected funnel; orthology confidence %s, no in-paralogue", orthology_confidence))]

## ---------------------------------------------------------------- release table
out <- m[, .(
  locus_tag, product, set, rank_published,
  # constraint
  piNpiS = round(piNpiS_shrunk, 4), constraint_pctile, callable_genomes = mean_cov_genomes,
  syn_seg, nonsyn_seg, segregating_sites = tot_seg, expected_variants = round(exp_tot, 1),
  # robustness
  pctile_median, pctile_q05, pctile_q95, rank_median, rank_q05, rank_q95,
  eligibility_retention = eligible_freq,
  # human homology
  homology_search_default = fifelse(is.na(hs_def_hit), "no hit", hs_def_hit),
  homology_search_sensitive = fifelse(is.na(hs_sen_hit), "no hit", hs_sen_hit),
  homology_sens_pident = hs_sen_pident, homology_sens_qcov = hs_sen_qcov,
  homology_sens_evalue = hs_sen_evalue,
  homology_threshold = "DIAMOND --very-sensitive, E<=1e-3 against the human reference proteome",
  homology_decision = fifelse(c4_no_human_sen, "pass: no detectable human homologue",
                              "cut: human homologue detectable"),
  # paralogy / orthology
  n_af_inparalogs, best_inparalog_acc, best_inparalog_pident,
  orthology_confidence, n_risk_flags, risk_flags, domain_mismatch_class,
  ct_orthologue = ct_locus_recomputed, ct_reciprocal,
  calb_orthologue = calb_locus_recomputed, calb_reciprocal,
  # structure
  mean_plddt, n_pockets, best_druggability, best_pocket_volume,
  best_pdb_id, best_pdb_identity = best_pdb, top_ligand,
  structure_caveat = "AlphaFold model and fpocket score are predicted druggability, not measured binding",
  # external validation
  calb_acc, calb_pident, calb_qcov, n_src_reporting, n_src_essential, sources_essential,
  external_validation_caveat = "published C. albicans orthologue essentiality; NOT an experimentally validated A. fumigatus target",
  ct_orthologue_crossgenus = ct_orthologue, ct_identity_pct, ct_piNpiS, ct_pctile,
  crossgenus_quadrant,
  # funnel status
  eligible_published, eligible_corrected, removal_reason,
  # triage
  triage, triage_rationale, confidence_rule)]
setorder(out, triage, -constraint_pctile)
fwrite(out, file.path(TAB, "rev_c2_candidate_evidence_master.csv"))

say("\n=== candidate evidence master: %d genes ===", nrow(out))
print(out[, .N, by = triage])
say("\n--- by triage class ---")
print(out[, .(locus_tag, product = substr(product, 1, 36), triage,
              hom = homology_decision, orth = orthology_confidence,
              inpara = n_af_inparalogs, ess = paste0(n_src_essential, "/", n_src_reporting))])
say("\nwrote rev_c2_candidate_evidence_master.csv")
