#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# Advisor round 2, WP1: one canonical filter table from which every published count derives.
#
# Why this exists. The funnel counts in the manuscript were hard-coded in
# analysis_code/ms_figures_tables.R, and the high-sensitivity human-homology re-screen was
# applied nowhere in code - the only artefact was a flat list of loci. The advisor's pass
# gate is that one reproducible table yields every published count and that no gene can be
# labelled host-filtered while retaining a disqualifying high-sensitivity hit.
#
# What it does. Recomputes the published funnel gene by gene from the canonical inputs,
# re-runs the WHOLE funnel with the sensitive homology screen substituted at the selectivity
# step - rather than subtracting the sensitive hits from the final 815, which is what the
# manuscript implicitly did - and writes a per-gene audit carrying, for every gene, its
# status under both rules, the first criterion it fails, and its human-homology evidence.
#
# Run: Rscript analysis_code/rev_c1_canonical_filter_audit.R
# ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table)})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
TD  <- "analysis/target_discovery"
TAB <- "output/table"

say <- function(...) cat(sprintf(...), "\n", sep = "")

## ---------------------------------------------------------------- canonical inputs
d <- fread(file.path(TD, "constraint_per_gene.tsv"))
stopifnot(nrow(d) == 9823, !any(duplicated(d$protein_id)))

# derived constraint statistics, identical to analysis_code/td_funnel.R
d[, tot_seg := syn_seg + nonsyn_seg]
d[, exp_nonsyn_frac := nonsyn_sites / (nonsyn_sites + syn_sites)]
d[, p_nonsyn_depl := ifelse(tot_seg >= 3, pbinom(nonsyn_seg, tot_seg, exp_nonsyn_frac), NA_real_)]
d[, piNpiS_shrunk := ((nonsyn_seg + 0.5) / nonsyn_sites) / ((syn_seg + 0.5) / syn_sites)]
cds_bp <- sum(d$nonsyn_sites + d$syn_sites)
dens   <- sum(d$tot_seg) / cds_bp
d[, exp_tot := (nonsyn_sites + syn_sites) * dens]
d[, p_tot_depl := ppois(tot_seg, exp_tot)]
# constraint percentile among all callable, variant-bearing genes (as in td_funnel.R)
ref_dist <- d[mean_cov_genomes >= 250 & tot_seg >= 3]$piNpiS_shrunk
d[, constraint_pctile := round(100 * sapply(piNpiS_shrunk,
        function(x) mean(ref_dist > x, na.rm = TRUE)), 1)]

HCOL <- c("protein_id", "human_hit", "pident", "len", "evalue", "bits", "qcov", "scov")
best_hit <- function(f) {
  h <- fread(f, col.names = HCOL)
  h[, .SD[which.min(evalue)], by = protein_id]
}
hum_def <- best_hit(file.path(TD, "hom/af293_vs_human.tsv"))       # DIAMOND default sensitivity
hum_sen <- best_hit(file.path(TD, "hom/af293_vs_human_sens.tsv"))  # DIAMOND --very-sensitive

par <- fread(file.path(TD, "hom/has_paralog.txt"), header = FALSE)$V1
xs  <- fread(file.path(TD, "hom/crossspecies_besthits.tsv"), header = FALSE,
             col.names = c("protein_id", "species"))
xsn <- xs[, .(n_species = uniqueN(species)), by = protein_id]

known <- c("AFUA_4G06890","AFUA_7G03740","AFUA_6G12400","AFUA_2G11010","AFUA_1G14870",
           "AFUA_2G03700","AFUA_1G14330","AFUA_2G02690","AFUA_5G04170","AFUA_5G09360",
           "AFUA_4G04820","AFUA_6G05300","AFUA_2G01260","AFUA_3G09850")

d <- merge(d, hum_def[, .(protein_id, hs_def_hit = human_hit, hs_def_pident = pident,
                          hs_def_qcov = qcov, hs_def_evalue = evalue, hs_def_bits = bits)],
           by = "protein_id", all.x = TRUE)
d <- merge(d, hum_sen[, .(protein_id, hs_sen_hit = human_hit, hs_sen_pident = pident,
                          hs_sen_qcov = qcov, hs_sen_evalue = evalue, hs_sen_bits = bits)],
           by = "protein_id", all.x = TRUE)
d <- merge(d, xsn, by = "protein_id", all.x = TRUE)
d[is.na(n_species), n_species := 0L]

## ---------------------------------------------------------------- per-gene criteria
# Each criterion is evaluated for every gene independently of the funnel order, so the
# table can answer "why was this gene dropped" and "what else would have dropped it".
d[, `:=`(
  c1_callable     = mean_cov_genomes >= 250,
  c2_evaluable    = tot_seg >= 3 | exp_tot >= 10,
  c3_constrained  = (!is.na(p_nonsyn_depl) & p_nonsyn_depl < 0.05) | p_tot_depl < 0.01,
  c4_no_human_def = is.na(hs_def_evalue),
  c4_no_human_sen = is.na(hs_sen_evalue),
  c5_no_paralog   = !(protein_id %in% par),
  c6_conserved    = n_species >= 8,
  c7_not_known    = !(locus_tag %in% known),
  c8_len_ok       = protein_len >= 100,
  c9_powered      = exp_tot >= 15
)]

# the sensitive screen must be a superset of the default screen; if it is not, the two
# searches are not comparable and nothing downstream can be trusted
viol <- d[c4_no_human_sen == TRUE & c4_no_human_def == FALSE]
if (nrow(viol)) {
  fwrite(viol[, .(protein_id, locus_tag, hs_def_hit, hs_def_evalue, hs_sen_hit, hs_sen_evalue)],
         file.path(TAB, "rev_c1_homology_superset_violations.csv"))
  say("WARNING: %d genes have a default-sensitivity hit but no sensitive hit", nrow(viol))
}

## ---------------------------------------------------------------- the funnel, both rules
STEPS <- c("c1_callable", "c2_evaluable", "c3_constrained", "HOM",
           "c5_no_paralog", "c6_conserved", "c7_not_known", "c8_len_ok", "c9_powered")
LABEL <- c("all Af293 protein-coding genes",
           "callable (aligned in >=250 of 315 genomes)",
           "evaluable (>=3 segregating sites or >=10 expected)",
           "constrained (nonsyn-depleted or variant-poor)",
           "no human homologue  <- selectivity",
           "no A. fumigatus paralogue  <- no redundancy",
           "conserved in >=8 of 9 other Aspergillus",
           "not an already-drugged target",
           "protein >=100 aa",
           "adequately powered (>=15 expected variants)")

run_funnel <- function(dt, hom_col) {
  keep <- rep(TRUE, nrow(dt)); counts <- nrow(dt)
  for (s in STEPS) {
    col  <- if (s == "HOM") hom_col else s
    keep <- keep & dt[[col]]
    counts <- c(counts, sum(keep))
  }
  list(keep = keep, counts = counts)
}
pub <- run_funnel(d, "c4_no_human_def")   # the published funnel
cor <- run_funnel(d, "c4_no_human_sen")   # the same funnel, sensitive screen substituted

d[, eligible_published := pub$keep]
d[, eligible_corrected := cor$keep]

# first criterion each gene fails, under the corrected rule
order_cor <- c("c1_callable", "c2_evaluable", "c3_constrained", "c4_no_human_sen",
               "c5_no_paralog", "c6_conserved", "c7_not_known", "c8_len_ok", "c9_powered")
names(LABEL) <- c("(none)", order_cor)
d[, removal_reason := {
  m <- as.matrix(d[, ..order_cor])
  idx <- apply(m, 1, function(r) { w <- which(!r); if (!length(w)) 0L else w[1] })
  ifelse(idx == 0L, "retained", order_cor[pmax(idx, 1)])
}]
d[eligible_corrected == TRUE, removal_reason := "retained"]

counts <- data.table(
  step              = seq_along(LABEL) - 1L,
  criterion         = c("(none)", order_cor),
  label             = unname(LABEL),
  n_published_rule  = pub$counts,
  n_corrected_rule  = cor$counts)
counts[, delta := n_corrected_rule - n_published_rule]
fwrite(counts, file.path(TAB, "rev_c1_funnel_counts_canonical.csv"))

say("\n=== FUNNEL: published (default-sensitivity homology) vs corrected (--very-sensitive) ===")
print(counts[, .(label, published = n_published_rule, corrected = n_corrected_rule, delta)])

n_pub_final <- tail(pub$counts, 1); n_cor_final <- tail(cor$counts, 1)

## ---------------------------------------------------------------- audit of the published 815
cand <- d[eligible_published == TRUE]
stopifnot(nrow(cand) == n_pub_final)
cand[, has_sensitive_human_hit := !c4_no_human_sen]
cand[, status_corrected := fifelse(eligible_corrected, "retained",
                            fifelse(!c4_no_human_sen, "cut: human homologue at high sensitivity",
                                    paste0("cut: ", removal_reason)))]
aud <- cand[, .(rank_published = NA_integer_, locus_tag, protein_id, product, protein_len,
                piNpiS_shrunk = round(piNpiS_shrunk, 4), exp_tot = round(exp_tot, 1),
                tot_seg, syn_seg, nonsyn_seg,
                human_hit_default = fifelse(is.na(hs_def_hit), "", hs_def_hit),
                human_hit_sensitive = fifelse(is.na(hs_sen_hit), "", hs_sen_hit),
                sens_pident = hs_sen_pident, sens_qcov = hs_sen_qcov,
                sens_evalue = hs_sen_evalue, sens_bits = hs_sen_bits,
                has_sensitive_human_hit, eligible_published, eligible_corrected,
                status_corrected)]
setorder(aud, -eligible_corrected, status_corrected, locus_tag)
fwrite(aud, file.path(TAB, "rev_c1_candidate_815_audit.csv"))

n_hit <- sum(cand$has_sensitive_human_hit)
say("\npublished candidate set                     %4d", nrow(cand))
say("  of these, carry a high-sensitivity hit     %4d (%.1f%%)", n_hit, 100 * n_hit / nrow(cand))
say("  survive the corrected funnel               %4d", sum(cand$eligible_corrected))
say("naive 'published minus sensitive hits'       %4d", nrow(cand) - n_hit)
say("corrected funnel run end to end              %4d", n_cor_final)
say("genes eligible ONLY under the corrected rule %4d",
    sum(d$eligible_corrected & !d$eligible_published))

## ---------------------------------------------------------------- Tier-1 re-evaluation
t1 <- fread(file.path(TAB, "ms_table3_tier1.csv"))
t1a <- merge(t1[, .(rank, locus_tag, product, plddt, pocket, pdb, pdb_identity, ligand,
                    human_homolog_sensitive_asreported = human_homolog_sensitive)],
             d[, .(locus_tag, piNpiS_shrunk, constraint_pctile, exp_tot, tot_seg,
                   syn_seg, nonsyn_seg, protein_len,
                   human_hit_sensitive = hs_sen_hit, sens_pident = hs_sen_pident,
                   sens_qcov = hs_sen_qcov, sens_evalue = hs_sen_evalue,
                   c4_no_human_sen, c5_no_paralog, c6_conserved, c9_powered,
                   eligible_published, eligible_corrected, removal_reason)],
             by = "locus_tag", all.x = TRUE)
setorder(t1a, rank)
fwrite(t1a, file.path(TAB, "rev_c1_tier1_reevaluation.csv"))
say("\nTier-1 candidates re-evaluated               %4d", nrow(t1a))
say("  retained under the corrected funnel        %4d", sum(t1a$eligible_corrected, na.rm = TRUE))
say("  carrying a high-sensitivity human hit      %4d", sum(!t1a$c4_no_human_sen, na.rm = TRUE))
if (any(!t1a$eligible_corrected, na.rm = TRUE))
  print(t1a[eligible_corrected == FALSE, .(rank, locus_tag, product, removal_reason)])

## ---------------------------------------------------------------- the full per-gene table
keep_cols <- c("protein_id", "locus_tag", "product", "protein_len", "n_codons",
               "syn_sites", "nonsyn_sites", "syn_seg", "nonsyn_seg", "stop_seg", "tot_seg",
               "piN_piS", "piNpiS_shrunk", "constraint_pctile", "exp_tot", "p_nonsyn_depl", "p_tot_depl",
               "mean_cov_genomes", "min_cov_genomes", "n_species",
               "hs_def_hit", "hs_def_pident", "hs_def_qcov", "hs_def_evalue",
               "hs_sen_hit", "hs_sen_pident", "hs_sen_qcov", "hs_sen_evalue",
               "c1_callable", "c2_evaluable", "c3_constrained", "c4_no_human_def",
               "c4_no_human_sen", "c5_no_paralog", "c6_conserved", "c7_not_known",
               "c8_len_ok", "c9_powered",
               "eligible_published", "eligible_corrected", "removal_reason")
fwrite(d[, ..keep_cols], file.path(TAB, "rev_c1_canonical_filter_audit.csv"))
say("\nwrote rev_c1_canonical_filter_audit.csv (%d genes), rev_c1_funnel_counts_canonical.csv,",
    nrow(d))
say("      rev_c1_candidate_815_audit.csv, rev_c1_tier1_reevaluation.csv")
