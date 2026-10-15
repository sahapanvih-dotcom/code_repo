#!/usr/bin/env Rscript
# =============================================================================
# rev_b3_pathway_enrichment.R
# Analysis B3 (manuscript revision): pathway-level functional enrichment of
# evolutionarily constrained genes, with an evaluable/callable/annotated
# background universe and matched null models.
#
# Why the tests are hand-rolled: topGO and clusterProfiler are NOT installed on
# this host and installation is not available, so Fisher/hypergeometric tests,
# covariate-adjusted logistic regression and matched-set permutation are
# implemented directly here. That is also what the advisor asked for, because
# an off-the-shelf tool would silently use "all genes in the genome" as the
# background, which is the wrong background for a constraint analysis: a gene
# can only be called constrained if it was callable AND evaluable.
#
# Annotation: KEGG REST (primary) + UniProtKB GO (secondary), both cached under
# analysis/enrichment/ on 2026-10-12. See analysis/enrichment/PROVENANCE.md.
#
# HONESTY RULE applied throughout: every tested term is reported, including the
# null results. Nothing is called "enriched" unless it survives BH FDR.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

B  <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
setwd(B)
ENR <- "analysis/enrichment"
# Advisor round 2 (WP5): the C. tropicalis estimator is a parameter. The project primary is
# the raw ratio on lineage representatives (piNpiS_dedup, see ctrop_constraint.py); the shrunk
# variant is a sensitivity analysis. The earlier run used the shrunk one without saying so.
CT_EST     <- Sys.getenv("REV_B3_CT_ESTIMATOR", unset = "piNpiS_dedup")
OUT_SUFFIX <- Sys.getenv("REV_B3_OUT_SUFFIX",   unset = "")
stopifnot(CT_EST %in% c("piNpiS_dedup", "piNpiS_shrunk_dedup"))
TAB <- "output/table"
FIG <- "output/figure"
RETRIEVAL_DATE <- "2026-10-12"   # annotation download date (see PROVENANCE.md)

NPERM      <- 10000L  # >= 1000 required; 10000 used (see note on the empirical floor)
MIN_TERM   <- 5L      # min genes of a term inside the universe to be testable
MAX_TERM   <- 1000L   # upper bound (drops the all-encompassing BRITE roots)
Q_STRONG   <- 0.10    # "strong constraint" = lowest decile of piN/piS
Q_SENS     <- 0.05    # sensitivity definition
Q_BOTH     <- 0.25    # for the both-genera intersection (decile^2 is too small)
FDR_ALPHA  <- 0.05
set.seed(20261012)

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# =============================================================================
# 1. ANNOTATION  (cached files only -- no network access at analysis time)
# =============================================================================
msg("== 1. loading cached annotation (retrieved %s) ==", RETRIEVAL_DATE)

read_link <- function(f, pfx) {
  x <- fread(file.path(ENR, f), header = FALSE, col.names = c("g", "t"))
  x[, .(gene = sub(paste0("^", pfx, ":"), "", g), term = sub("^(path|md|br|ko):", "", t))]
}

kegg_path_afm <- read_link("kegg_afm_gene_pathway.tsv", "afm")
kegg_mod_afm  <- read_link("kegg_afm_module.tsv",       "afm")
kegg_bri_afm  <- read_link("kegg_afm_brite.tsv",        "afm")
kegg_path_ctp <- read_link("kegg_ctp_gene_pathway.tsv", "ctp")
kegg_mod_ctp  <- read_link("kegg_ctp_module.tsv",       "ctp")

# ---- term names -------------------------------------------------------------
pname_afm <- fread(file.path(ENR, "kegg_afm_pathway_names.tsv"), header = FALSE,
                   col.names = c("term", "name"))
pname_ctp <- fread(file.path(ENR, "kegg_ctp_pathway_names.tsv"), header = FALSE,
                   col.names = c("term", "name"))
pname_afm[, name := sub(" - Aspergillus fumigatus$", "", name)]
pname_ctp[, name := sub(" - Candida tropicalis$",   "", name)]
mname <- fread(file.path(ENR, "kegg_module_names.tsv"), header = FALSE,
               col.names = c("term", "name"))
mname[, term := paste0("afm_", term)]
mname2 <- copy(mname); mname2[, term := sub("^afm_", "ctp_", term)]
mod_names <- rbind(mname, mname2)
bname <- fread(file.path(ENR, "kegg_brite_names.tsv"), header = FALSE,
               col.names = c("term", "name"))
bname[, term := sub("^ko|^br", "afm", term)]

# ---- GO (reachable on 2026-10-12; see PROVENANCE.md) ------------------------
go_afm_raw <- fread(file.path(ENR, "uniprot_afm_go.tsv"), sep = "\t", quote = "")
setnames(go_afm_raw, c("uniprot", "oln", "gene_primary", "go_ids", "go_full"))
go_ctp_raw <- fread(file.path(ENR, "uniprot_ctp_go.tsv"), sep = "\t", quote = "")
setnames(go_ctp_raw, c("uniprot", "oln", "gene_names", "go_ids", "go_full"))

afua2up <- fread("analysis/target_discovery/struct/afua_to_uniprot.tsv")
go_afm <- merge(afua2up, go_afm_raw[go_ids != "", .(uniprot, go_ids)], by = "uniprot")
go_afm <- go_afm[, .(term = trimws(unlist(strsplit(go_ids, ";")))), by = .(gene = locus_tag)]
go_afm <- unique(go_afm[term != ""])

go_ctp_raw[, ctrg := {
  m <- regmatches(gene_names, regexpr("CTRG_[0-9]+", gene_names)); ifelse(lengths(m) == 0, NA, m)
}]
go_ctp <- go_ctp_raw[go_ids != "" & !is.na(ctrg)]
go_ctp <- go_ctp[, .(term = trimws(unlist(strsplit(go_ids, ";")))), by = .(gene = ctrg)]
go_ctp <- unique(go_ctp[term != ""])

# GO id -> name, harvested from the same cached UniProt files (no term invented)
go_names <- rbind(
  go_afm_raw[go_full != "", .(x = unlist(strsplit(go_full, "; ")))],
  go_ctp_raw[go_full != "", .(x = unlist(strsplit(go_full, "; ")))])
go_names <- unique(go_names[grepl("\\[GO:[0-9]+\\]$", x)])
go_names[, term := sub(".*\\[(GO:[0-9]+)\\]$", "\\1", x)]
go_names[, name := sub(" \\[GO:[0-9]+\\]$", "", x)]
go_names <- unique(go_names[, .(term, name)])[!duplicated(term)]

term_names <- unique(rbind(pname_afm, pname_ctp, mod_names, bname, go_names))
term_names <- term_names[!duplicated(term)]

msg("   KEGG pathway links: afm=%d genes / %d pathways ; ctp=%d genes / %d pathways",
    uniqueN(kegg_path_afm$gene), uniqueN(kegg_path_afm$term),
    uniqueN(kegg_path_ctp$gene), uniqueN(kegg_path_ctp$term))
msg("   GO links: afm=%d genes / %d terms ; ctp=%d genes / %d terms",
    uniqueN(go_afm$gene), uniqueN(go_afm$term), uniqueN(go_ctp$gene), uniqueN(go_ctp$term))

# =============================================================================
# 2. CONSTRAINT TABLES AND THE GENE UNIVERSE
#    Universe = EVALUABLE & CALLABLE & ANNOTATED, never "all genes in genome".
# =============================================================================
msg("== 2. building gene universe ==")
univ_steps <- list()
add_step <- function(sp, label, n) univ_steps[[length(univ_steps) + 1L]] <<-
  data.table(species = sp, step = label, n_genes = n)

# ---------------- A. fumigatus ----------------------------------------------
af <- fread("analysis/target_discovery/constraint_per_gene_dedup.tsv")
add_step("A_fumigatus", "1. all Af293 protein-coding genes in constraint table", nrow(af))
af[, tot_seg := syn_seg + nonsyn_seg]
af[, cds_bp  := syn_sites + nonsyn_sites]
dens_af <- sum(af$tot_seg) / sum(af$cds_bp)
af[, exp_tot := cds_bp * dens_af]
# shrunk piN/piS: defined even when syn_seg is small (same estimator as td_funnel.R)
af[, piNpiS_shrunk := ((nonsyn_seg + 0.5) / nonsyn_sites) / ((syn_seg + 0.5) / syn_sites)]

# callable: aligned in >= 80% of the 187 de-duplicated genomes (cov >= 150).
# The published funnel used >=250 of 315 raw genomes (79.4%); 150/187 = 80.2% is
# the same criterion transferred to the de-duplicated panel used here.
AF_COV_MIN <- 150
af_c <- af[mean_cov_genomes >= AF_COV_MIN]
add_step("A_fumigatus", sprintf("2. callable: aligned in >=%d of 187 dedup genomes", AF_COV_MIN), nrow(af_c))
# evaluable: enough observed variants, or enough expected variants to detect depletion
af_e <- af_c[tot_seg >= 3 | exp_tot >= 10]
add_step("A_fumigatus", "3. evaluable: tot_seg>=3 or expected variants>=10", nrow(af_e))
af_e <- af_e[is.finite(piNpiS_shrunk) & protein_len > 0]
add_step("A_fumigatus", "4. finite piN/piS and protein_len>0 (constraint-evaluable universe)", nrow(af_e))

# covariates
xs <- fread("analysis/target_discovery/hom/crossspecies_besthits.tsv", header = FALSE,
            col.names = c("protein_id", "species"))
xsn <- xs[, .(n_ortho = uniqueN(species)), by = protein_id]
af_e <- merge(af_e, xsn, by = "protein_id", all.x = TRUE)
af_e[is.na(n_ortho), n_ortho := 0L]
hum_sens <- fread("analysis/target_discovery/hom/human_homolog_sensitive_loci.txt",
                  header = FALSE)$V1
hum_strict <- fread("analysis/target_discovery/hom/af293_vs_human.tsv", header = FALSE,
                    col.names = c("protein_id", "hit", "pident", "alen", "evalue",
                                  "bits", "qcov", "scov"))
hum_strict_best <- hum_strict[, .SD[which.min(evalue)], by = protein_id]
af_e[, has_human_homolog_sens   := locus_tag %in% hum_sens]
af_e <- merge(af_e, hum_strict_best[, .(protein_id, human_pident = pident,
                                        human_qcov = qcov, human_evalue = evalue)],
              by = "protein_id", all.x = TRUE)
af_e[, has_human_homolog_strict := !is.na(human_evalue)]

# ---------------- C. tropicalis ---------------------------------------------
ct <- fread("analysis/candida/constraint/ctrop_constraint_per_gene.tsv")
add_step("C_tropicalis", "1. all MYA-3404 protein-coding genes in constraint table", nrow(ct))
ct[, tot_seg := syn_seg + nonsyn_seg]
ct[, cds_bp  := syn_sites + nonsyn_sites]
dens_ct <- sum(ct$tot_seg) / sum(ct$cds_bp)
ct[, exp_tot := cds_bp * dens_ct]
CT_COV_MIN <- 340                       # 340 / 424 isolates = 80.2%, same rule as above
ct_c <- ct[n_iso_callable >= CT_COV_MIN]
add_step("C_tropicalis", sprintf("2. callable: >=%d of 424 isolates callable", CT_COV_MIN), nrow(ct_c))
ct_e <- ct_c[tot_seg >= 3 | exp_tot >= 10]
add_step("C_tropicalis", "3. evaluable: tot_seg>=3 or expected variants>=10", nrow(ct_e))
ct_e <- ct_e[is.finite(get(CT_EST)) & protein_len > 0]
add_step("C_tropicalis", "4. finite piN/piS and protein_len>0 (constraint-evaluable universe)", nrow(ct_e))

# number of detectable cross-genus orthologues for a C. tropicalis gene (0-2:
# A. fumigatus and/or C. albicans), from the cross-genus ortholog map
om <- fread("analysis/crossgenus/ortho_map.tsv", na.strings = c("", "NA"))
ct_orth <- om[!is.na(ct_locus), .(n_ortho = as.integer(1L + any(!is.na(calb_acc)))),
              by = .(locus_tag = ct_locus)]
ct_e <- merge(ct_e, ct_orth, by = "locus_tag", all.x = TRUE)
ct_e[is.na(n_ortho), n_ortho := 0L]

# ---------------- cross-genus ortholog-pair universe ------------------------
pairs <- om[!is.na(ct_locus), .(af_locus, ct_locus, ct_pident)]
pairs <- pairs[af_locus %in% af_e$locus_tag & ct_locus %in% ct_e$locus_tag]
add_step("both_genera", "1. A.fum-C.trop ortholog pairs, both genes constraint-evaluable",
         nrow(pairs))

# ---------------- annotated universes ---------------------------------------
ann_sets <- list(
  list(set = "A_fumigatus", source = "KEGG_PATHWAY", t2g = kegg_path_afm),
  list(set = "A_fumigatus", source = "KEGG_MODULE",  t2g = kegg_mod_afm),
  list(set = "A_fumigatus", source = "KEGG_BRITE",   t2g = kegg_bri_afm),
  list(set = "A_fumigatus", source = "GO",           t2g = go_afm),
  list(set = "C_tropicalis", source = "KEGG_PATHWAY", t2g = kegg_path_ctp),
  list(set = "C_tropicalis", source = "KEGG_MODULE",  t2g = kegg_mod_ctp),
  list(set = "C_tropicalis", source = "GO",           t2g = go_ctp))
for (a in ann_sets) {
  pool <- if (a$set == "A_fumigatus") af_e$locus_tag else ct_e$locus_tag
  add_step(a$set, sprintf("5. + annotated in %s (universe used for %s tests)", a$source, a$source),
           sum(pool %in% a$t2g$gene))
}
add_step("both_genera", "2. + A. fumigatus gene annotated in KEGG_PATHWAY",
         sum(pairs$af_locus %in% kegg_path_afm$gene))
add_step("both_genera", "3. + A. fumigatus gene annotated in GO",
         sum(pairs$af_locus %in% go_afm$gene))

univ_tab <- rbindlist(univ_steps)
msg("   A. fumigatus constraint-evaluable universe: %d genes", nrow(af_e))
msg("   C. tropicalis constraint-evaluable universe: %d genes", nrow(ct_e))
msg("   cross-genus ortholog-pair universe: %d pairs", nrow(pairs))

# =============================================================================
# 3. GENE SETS
# =============================================================================
msg("== 3. defining constrained gene sets ==")
af_thr10 <- quantile(af_e$piNpiS_shrunk, Q_STRONG, na.rm = TRUE)
af_thr05 <- quantile(af_e$piNpiS_shrunk, Q_SENS,   na.rm = TRUE)
af_thr25 <- quantile(af_e$piNpiS_shrunk, Q_BOTH,   na.rm = TRUE)
ct_thr10 <- quantile(ct_e[[CT_EST]], Q_STRONG, na.rm = TRUE)
ct_thr05 <- quantile(ct_e[[CT_EST]], Q_SENS,   na.rm = TRUE)
ct_thr25 <- quantile(ct_e[[CT_EST]], Q_BOTH,   na.rm = TRUE)

af_e[, `:=`(constr10 = piNpiS_shrunk <= af_thr10,
            constr05 = piNpiS_shrunk <= af_thr05,
            constr25 = piNpiS_shrunk <= af_thr25)]
ct_e[, `:=`(constr10 = get(CT_EST) <= ct_thr10,
            constr05 = get(CT_EST) <= ct_thr05,
            constr25 = get(CT_EST) <= ct_thr25)]

pairs <- merge(pairs, af_e[, .(af_locus = locus_tag, af_c10 = constr10, af_c25 = constr25,
                               af_len = cds_bp, af_n_ortho = n_ortho)], by = "af_locus")
pairs <- merge(pairs, ct_e[, .(ct_locus = locus_tag, ct_c10 = constr10, ct_c25 = constr25)],
               by = "ct_locus")
pairs[, both10 := af_c10 & ct_c10]
pairs[, both25 := af_c25 & ct_c25]

gene_sets_tab <- data.table(
  gene_set = c("AF_constrained_q10", "AF_constrained_q05",
               "CT_constrained_q10", "CT_constrained_q05",
               "BOTH_constrained_q25", "BOTH_constrained_q10",
               "NEGCTRL_AF_matched_random", "NEGCTRL_AF_uniform_random"),
  definition = c(
    sprintf("A. fumigatus piNpiS_shrunk <= %.4f (lowest %.0f%% of the evaluable universe)", af_thr10, 100*Q_STRONG),
    sprintf("A. fumigatus piNpiS_shrunk <= %.4f (lowest %.0f%%)", af_thr05, 100*Q_SENS),
    sprintf("C. tropicalis %s <= %.4f (lowest %.0f%%)", CT_EST, ct_thr10, 100*Q_STRONG),
    sprintf("C. tropicalis %s <= %.4f (lowest %.0f%%)", CT_EST, ct_thr05, 100*Q_SENS),
    sprintf("ortholog pairs in the lowest %.0f%% of BOTH genera", 100*Q_BOTH),
    "ortholog pairs in the lowest 10% of BOTH genera (strict; reported for completeness)",
    "negative control 1: one draw from the same matched null used by the permutation test (strata = length decile x annotation-depth bin); calibrates the permutation framework",
    "negative control 2: uniform random draw from the annotated universe; calibrates the unadjusted Fisher test"),
  n_genes = c(sum(af_e$constr10), sum(af_e$constr05), sum(ct_e$constr10), sum(ct_e$constr05),
              sum(pairs$both25), sum(pairs$both10), sum(af_e$constr10), sum(af_e$constr10)))
print(gene_sets_tab[, .(gene_set, n_genes)])
msg("   expected both-genera overlap under independence: q25 %.0f, q10 %.0f",
    nrow(pairs) * Q_BOTH^2, nrow(pairs) * Q_STRONG^2)

# excess overlap of constraint across genera (a prerequisite for the BOTH tests)
ov_tab <- rbindlist(lapply(c("10", "25"), function(q) {
  a <- pairs[[paste0("af_c", q)]]; b <- pairs[[paste0("ct_c", q)]]
  ft <- fisher.test(table(factor(a, c(FALSE, TRUE)), factor(b, c(FALSE, TRUE))))
  data.table(quantile_cut = paste0("q", q), n_pairs = length(a),
             n_af = sum(a), n_ct = sum(b), n_both = sum(a & b),
             expected_both = length(a) * mean(a) * mean(b),
             odds_ratio = unname(ft$estimate), ci_lo = ft$conf.int[1],
             ci_hi = ft$conf.int[2], p_value = ft$p.value)
}))

# =============================================================================
# 4. TEST ENGINE
# =============================================================================
# (a) Fisher exact / hypergeometric on the annotated universe
# (b) logistic regression: constrained ~ in_term + log10(CDS bp)
#                                      + annotation depth + n detectable orthologues
# (c) matched-set permutation: NPERM null "constrained" sets of identical size,
#     drawn within strata of (log10 CDS-length decile) x (annotation-depth bin).
#     Matching tolerance = same length decile (so |delta log10 len| is bounded by
#     the width of one decile of the universe) and identical annotation-depth bin.
# =============================================================================

make_strata <- function(len_bp, ann_depth) {
  lg <- log10(len_bp)
  qs <- unique(quantile(lg, probs = seq(0, 1, 0.1), na.rm = TRUE))
  lbin <- cut(lg, breaks = qs, include.lowest = TRUE, labels = FALSE)
  abin <- ifelse(ann_depth <= 1, 1L, ifelse(ann_depth <= 3, 2L, 3L))
  list(stratum = paste0("L", lbin, "_A", abin), lbin = lbin, abin = abin,
       decile_width = median(diff(qs)))
}

# Returns a closure that draws one matched null set: the same number of genes per
# stratum as the observed constrained set, sampled without replacement from the
# whole universe within that stratum. Strata are precomputed once for speed.
make_sampler <- function(stratum, hit) {
  idx_by <- split(seq_along(stratum), stratum)
  need   <- table(stratum[hit])
  pools  <- idx_by[names(need)]
  ks     <- as.integer(need)
  function() {
    unlist(lapply(seq_along(pools), function(i) {
      p <- pools[[i]]; k <- ks[i]
      if (length(p) <= k) p else p[sample.int(length(p), k)]
    }), use.names = FALSE)
  }
}

run_tests <- function(label, source_name, genes, hit, len_bp, n_ortho, t2g, nperm = NPERM) {
  t2g <- unique(t2g[gene %in% genes, .(gene, term)])
  ann_depth <- t2g[, .N, by = gene][match(genes, gene)]$N
  ann_depth[is.na(ann_depth)] <- 0L
  keep <- ann_depth > 0                       # ANNOTATED universe
  genes <- genes[keep]; hit <- hit[keep]; len_bp <- len_bp[keep]
  n_ortho <- n_ortho[keep]; ann_depth <- ann_depth[keep]
  t2g <- t2g[gene %in% genes]
  N <- length(genes); n_set <- sum(hit)

  terms <- t2g[, .N, by = term][N >= MIN_TERM & N <= MAX_TERM]$term
  t2g <- t2g[term %in% terms]
  if (!length(terms)) return(NULL)
  M <- sparseMatrix(i = match(t2g$gene, genes), j = match(t2g$term, terms),
                    x = 1, dims = c(N, length(terms)), dimnames = list(genes, terms))

  n_term <- Matrix::colSums(M)
  obs    <- as.numeric(Matrix::crossprod(M, as.numeric(hit)))

  # ---- (a) Fisher exact ----------------------------------------------------
  fi <- rbindlist(lapply(seq_along(terms), function(j) {
    a <- obs[j]; b <- n_set - a; c <- n_term[j] - a; d <- N - n_set - c
    ft <- fisher.test(matrix(c(a, c, b, d), 2, 2))
    ha <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))   # Haldane-Anscombe OR
    data.table(or_fisher = unname(ft$estimate), or_lo = ft$conf.int[1],
               or_hi = ft$conf.int[2], or_haldane = ha, p_fisher = ft$p.value)
  }))

  # ---- (b) logistic regression with confounder covariates -------------------
  loglen <- log10(len_bp); lad <- log10(ann_depth)
  lo <- rbindlist(lapply(seq_along(terms), function(j) {
    inT <- as.numeric(M[, j])
    if (obs[j] == 0 || obs[j] == n_term[j] || obs[j] == n_set)
      return(data.table(beta = NA_real_, beta_lo = NA_real_, beta_hi = NA_real_,
                        p_logit = NA_real_, logit_note = "separable (0 or complete overlap)"))
    fit <- try(suppressWarnings(glm(hit ~ inT + loglen + lad + n_ortho,
                                    family = binomial())), silent = TRUE)
    if (inherits(fit, "try-error"))
      return(data.table(beta = NA_real_, beta_lo = NA_real_, beta_hi = NA_real_,
                        p_logit = NA_real_, logit_note = "glm failed"))
    s <- summary(fit)$coefficients
    if (!"inT" %in% rownames(s))
      return(data.table(beta = NA_real_, beta_lo = NA_real_, beta_hi = NA_real_,
                        p_logit = NA_real_, logit_note = "term dropped"))
    b <- s["inT", 1]; se <- s["inT", 2]
    data.table(beta = b, beta_lo = b - 1.96 * se, beta_hi = b + 1.96 * se,
               p_logit = s["inT", 4],
               logit_note = if (!fit$converged) "glm did not converge" else "ok")
  }))

  # ---- (c) matched-set permutation ----------------------------------------
  st <- make_strata(len_bp, ann_depth)
  draw <- make_sampler(st$stratum, hit)
  nT <- length(terms)
  cnt_ge <- integer(nT); cnt_le <- integer(nT)
  s1 <- numeric(nT); s2 <- numeric(nT)
  obs_ll <- mean(loglen[hit]); obs_ad <- mean(ann_depth[hit])
  sum_dl <- 0; sum_da <- 0; max_dl <- 0; max_da <- 0
  for (p in seq_len(nperm)) {
    idx <- draw()
    v <- numeric(N); v[idx] <- 1
    x <- as.numeric(Matrix::crossprod(M, v))
    cnt_ge <- cnt_ge + (x >= obs); cnt_le <- cnt_le + (x <= obs)
    s1 <- s1 + x; s2 <- s2 + x * x
    dl <- mean(loglen[idx]) - obs_ll; da <- mean(ann_depth[idx]) - obs_ad
    sum_dl <- sum_dl + dl; sum_da <- sum_da + da
    max_dl <- max(max_dl, abs(dl)); max_da <- max(max_da, abs(da))
  }
  nm <- s1 / nperm
  ns <- sqrt(pmax(s2 / nperm - nm^2, 0) * nperm / (nperm - 1))
  pe <- (1 + cnt_ge) / (nperm + 1)
  pd <- (1 + cnt_le) / (nperm + 1)
  pm <- pmin(1, 2 * pmin(pe, pd))
  # The empirical permutation p-value cannot go below 2/(nperm+1). In a family of
  # ~1000 terms that floor is above the BH threshold for the top-ranked term, so a
  # purely empirical permutation p can never reach FDR significance however real
  # the effect is. We therefore ALSO report a z-score p-value computed against the
  # same matched null (mean and sd of the matched null count); that statistic has
  # no resolution floor and is what the permutation-based FDR uses.
  pz <- 2 * pnorm(-abs((obs - nm) / ifelse(ns > 0, ns, NA_real_)))

  res <- data.table(
    gene_set = label, annotation_source = source_name, term = terms,
    term_name = term_names$name[match(terms, term_names$term)],
    n_universe = N, n_universe_in_term = as.integer(n_term),
    n_set = n_set, n_set_in_term = as.integer(obs),
    expected_in_term_hypergeom = n_set * n_term / N,
    fold_enrichment = (obs / n_set) / (n_term / N))
  res <- cbind(res, fi, lo)
  res[, `:=`(perm_null_mean = nm, perm_null_sd = ns,
             perm_fold = obs / pmax(nm, 1e-9),
             perm_z = (obs - nm) / ifelse(ns > 0, ns, NA_real_),
             p_perm_enrich = pe, p_perm_deplete = pd, p_perm_two_sided = pm,
             p_perm_z = pz, n_perm = nperm,
             perm_empirical_floor = 2 / (nperm + 1),
             perm_at_empirical_floor = pm <= 2 / (nperm + 1) + 1e-12)]
  attr(res, "matching") <- data.table(
    gene_set = label, annotation_source = source_name, n_strata = uniqueN(st$stratum),
    n_perm = nperm, length_decile_width_log10 = st$decile_width,
    mean_delta_log10_len_null_minus_obs = sum_dl / nperm,
    max_abs_delta_log10_len = max_dl,
    mean_delta_annotation_depth = sum_da / nperm,
    max_abs_delta_annotation_depth = max_da)
  res
}

# ---- assemble the jobs ------------------------------------------------------
msg("== 4. running tests (Fisher + logistic + %d matched permutations) ==", NPERM)

# Negative controls, both of the same size as AF_constrained_q10:
#  (1) one single draw from the SAME matched null used by the permutation test
#      (strata = length decile x KEGG annotation-depth bin). This is the control
#      that tells us whether the permutation framework is calibrated.
#  (2) a plain uniform random draw from the annotated universe. This is the
#      control that tells us whether the unadjusted Fisher test is calibrated.
# Both controls are drawn over the FULL A. fumigatus constraint-evaluable universe
# (n = 9067) so that they have the same size as AF_constrained_q10 and -- because
# annotation status is part of the matching -- the same per-annotation-source
# subset sizes. Matching strata: log10 CDS-length decile x annotation-status cell
# (annotated in KEGG pathways yes/no) x (annotated in GO yes/no).
af_e[, ann_kegg := locus_tag %in% kegg_path_afm$gene]
af_e[, ann_go   := locus_tag %in% go_afm$gene]
lgq <- unique(quantile(log10(af_e$cds_bp), seq(0, 1, 0.1)))
neg_stratum <- paste0("L", cut(log10(af_e$cds_bp), lgq, include.lowest = TRUE, labels = FALSE),
                      "_K", as.integer(af_e$ann_kegg), "_G", as.integer(af_e$ann_go))
draw_neg <- make_sampler(neg_stratum, af_e$constr10)
neg_idx  <- draw_neg()
af_e[, negctrl := FALSE][neg_idx, negctrl := TRUE]
n_target <- sum(af_e$constr10)
af_e[, negctrl_unif := FALSE][sample.int(nrow(af_e), n_target), negctrl_unif := TRUE]
msg("   negative controls: matched draw %d genes (%.0f%% overlap with the real set), uniform draw %d genes",
    sum(af_e$negctrl), 100 * mean(af_e$negctrl[af_e$constr10]), sum(af_e$negctrl_unif))

jobs <- list()
for (src in c("KEGG_PATHWAY", "KEGG_MODULE", "KEGG_BRITE", "GO")) {
  t2g <- switch(src, KEGG_PATHWAY = kegg_path_afm, KEGG_MODULE = kegg_mod_afm,
                KEGG_BRITE = kegg_bri_afm, GO = go_afm)
  for (s in c("AF_constrained_q10", "AF_constrained_q05",
              "NEGCTRL_AF_matched_random", "NEGCTRL_AF_uniform_random")) {
    hitv <- switch(s, AF_constrained_q10 = af_e$constr10,
                   AF_constrained_q05 = af_e$constr05,
                   NEGCTRL_AF_matched_random = af_e$negctrl, af_e$negctrl_unif)
    jobs[[length(jobs) + 1L]] <- list(label = s, source = src, genes = af_e$locus_tag,
                                      hit = hitv, len = af_e$cds_bp,
                                      n_ortho = af_e$n_ortho, t2g = t2g)
  }
}
for (src in c("KEGG_PATHWAY", "KEGG_MODULE", "GO")) {
  t2g <- switch(src, KEGG_PATHWAY = kegg_path_ctp, KEGG_MODULE = kegg_mod_ctp, GO = go_ctp)
  for (s in c("CT_constrained_q10", "CT_constrained_q05")) {
    hitv <- if (s == "CT_constrained_q10") ct_e$constr10 else ct_e$constr05
    jobs[[length(jobs) + 1L]] <- list(label = s, source = src, genes = ct_e$locus_tag,
                                      hit = hitv, len = ct_e$cds_bp,
                                      n_ortho = ct_e$n_ortho, t2g = t2g)
  }
}
for (src in c("KEGG_PATHWAY", "KEGG_MODULE", "GO")) {
  t2g <- switch(src, KEGG_PATHWAY = kegg_path_afm, KEGG_MODULE = kegg_mod_afm, GO = go_afm)
  for (s in c("BOTH_constrained_q25", "BOTH_constrained_q10")) {
    hitv <- if (s == "BOTH_constrained_q25") pairs$both25 else pairs$both10
    jobs[[length(jobs) + 1L]] <- list(label = s, source = src, genes = pairs$af_locus,
                                      hit = hitv, len = pairs$af_len,
                                      n_ortho = pairs$af_n_ortho, t2g = t2g)
  }
}

out <- list(); matchq <- list()
for (j in jobs) {
  r <- run_tests(j$label, j$source, j$genes, j$hit, j$len, j$n_ortho, j$t2g)
  if (is.null(r)) { msg("   %-26s %-13s no testable terms", j$label, j$source); next }
  msg("   %-26s %-13s %4d terms | universe %5d | set %4d",
      j$label, j$source, nrow(r), r$n_universe[1], r$n_set[1])
  matchq[[length(matchq) + 1L]] <- attr(r, "matching")
  out[[length(out) + 1L]] <- r
}
res <- rbindlist(out)
match_tab <- rbindlist(matchq)

# ---- replicate calibration: many independent negative-control draws ---------
# A single negative-control set only gives one realisation. To state whether the
# tests are CALIBRATED we repeat the control draw NCALIB times and record how
# many terms each replicate declares significant after BH. A calibrated test
# should declare ~0 in nearly every replicate.
NCALIB <- 20L
fisher_only <- function(genes, hit, t2g) {
  t2g <- unique(t2g[gene %in% genes, .(gene, term)])
  ad <- t2g[, .N, by = gene][match(genes, gene)]$N; ad[is.na(ad)] <- 0L
  k <- ad > 0; genes <- genes[k]; hit <- hit[k]
  t2g <- t2g[gene %in% genes]
  N <- length(genes); n_set <- sum(hit)
  tm <- t2g[, .N, by = term][N >= MIN_TERM & N <= MAX_TERM]$term
  t2g <- t2g[term %in% tm]
  M <- sparseMatrix(i = match(t2g$gene, genes), j = match(t2g$term, tm),
                    x = 1, dims = c(N, length(tm)))
  nt <- Matrix::colSums(M); ob <- as.numeric(Matrix::crossprod(M, as.numeric(hit)))
  p <- vapply(seq_along(tm), function(j) {
    a <- ob[j]; b <- n_set - a; c <- nt[j] - a; d <- N - n_set - c
    fisher.test(matrix(c(a, c, b, d), 2, 2))$p.value }, numeric(1))
  list(p = p, n_terms = length(tm), n_set_annotated = n_set)
}
calib_rows <- list()
for (src in c("KEGG_PATHWAY", "GO")) {
  t2g <- if (src == "KEGG_PATHWAY") kegg_path_afm else go_afm
  for (mode in c("uniform", "matched")) {
    for (r in seq_len(NCALIB)) {
      hv <- rep(FALSE, nrow(af_e))
      if (mode == "uniform") hv[sample.int(nrow(af_e), n_target)] <- TRUE
      else hv[draw_neg()] <- TRUE
      fo <- fisher_only(af_e$locus_tag, hv, t2g)
      q <- p.adjust(fo$p, "BH")
      calib_rows[[length(calib_rows) + 1L]] <- data.table(
        annotation_source = src, control_mode = mode, replicate = r,
        n_terms = fo$n_terms, n_set_annotated = fo$n_set_annotated,
        min_p_fisher = min(fo$p), n_BH_significant = sum(q < FDR_ALPHA),
        frac_p_below_0.05 = mean(fo$p < 0.05))
    }
  }
}
calib_tab <- rbindlist(calib_rows)
calib_sum <- calib_tab[, .(n_replicates = .N, n_terms = n_terms[1],
                           median_n_BH_sig = median(n_BH_significant),
                           mean_n_BH_sig = mean(n_BH_significant),
                           max_n_BH_sig = max(n_BH_significant),
                           frac_replicates_with_any_hit = mean(n_BH_significant > 0),
                           median_frac_p_below_0.05 = median(frac_p_below_0.05)),
                       by = .(annotation_source, control_mode)]
msg("== 4b. replicate calibration (%d draws per cell) ==", NCALIB)
print(calib_sum)

# =============================================================================
# 5. BH FDR
#    Family = gene_set x annotation_source x method (the natural testing family:
#    one annotation system interrogated with one statistic for one gene set).
#    A global BH across every test in the study is also reported.
# =============================================================================
msg("== 5. Benjamini-Hochberg FDR ==")
res[, q_fisher := p.adjust(p_fisher, "BH"), by = .(gene_set, annotation_source)]
res[, q_logit  := p.adjust(p_logit,  "BH"), by = .(gene_set, annotation_source)]
res[, q_perm   := p.adjust(p_perm_two_sided, "BH"), by = .(gene_set, annotation_source)]
res[, q_perm_z := p.adjust(p_perm_z, "BH"), by = .(gene_set, annotation_source)]
res[, q_fisher_global := p.adjust(p_fisher, "BH")]
res[, q_logit_global  := p.adjust(p_logit,  "BH")]
res[, q_perm_global   := p.adjust(p_perm_two_sided, "BH")]
res[, q_perm_z_global := p.adjust(p_perm_z, "BH")]
res[, direction := ifelse(n_set_in_term > expected_in_term_hypergeom, "enriched", "depleted")]
# NOTE ON WHICH PERMUTATION P-VALUE COUNTS.
# The empirical matched-permutation p (p_perm_two_sided) is the headline statistic:
# it is exact-by-construction and both negative controls give zero BH-significant
# terms with it. Its only drawback is a resolution floor of 2/(nperm+1) = 2e-4.
# The z-score version (p_perm_z) has no floor but is ANTI-CONSERVATIVE for sparse
# terms -- in the negative controls it produces false positives among GO terms
# whose matched-null mean count is < 2. It is therefore kept in the output as a
# ranking diagnostic only and is excluded from the significance flags below.
res[, significant_any := (!is.na(q_fisher) & q_fisher < FDR_ALPHA) |
      (!is.na(q_logit) & q_logit < FDR_ALPHA) |
      (!is.na(q_perm) & q_perm < FDR_ALPHA)]
res[, significant_all_three := (!is.na(q_fisher) & q_fisher < FDR_ALPHA) &
      (!is.na(q_logit) & q_logit < FDR_ALPHA) &
      (!is.na(q_perm) & q_perm < FDR_ALPHA)]
res[, perm_z_is_diagnostic_only := TRUE]

sig_counts <- res[, .(n_terms_tested = .N,
                      n_q_fisher_sig = sum(q_fisher < FDR_ALPHA, na.rm = TRUE),
                      n_q_logit_sig  = sum(q_logit  < FDR_ALPHA, na.rm = TRUE),
                      n_q_perm_empirical_sig = sum(q_perm < FDR_ALPHA, na.rm = TRUE),
                      n_q_perm_z_sig = sum(q_perm_z < FDR_ALPHA, na.rm = TRUE),
                      n_sig_all_three_methods = sum(significant_all_three, na.rm = TRUE),
                      n_q_fisher_sig_global = sum(q_fisher_global < FDR_ALPHA, na.rm = TRUE),
                      min_p_fisher = min(p_fisher, na.rm = TRUE)),
                  by = .(gene_set, annotation_source)]
print(sig_counts)

# agreement between the covariate-adjusted model and the matched permutation
agree <- res[!is.na(p_logit) & !is.na(p_perm_z),
             .(n = .N,
               spearman_logit_vs_perm_z = suppressWarnings(
                 cor(-log10(p_logit), -log10(p_perm_z), method = "spearman")),
               spearman_fisher_vs_perm_z = suppressWarnings(
                 cor(-log10(p_fisher), -log10(p_perm_z), method = "spearman")),
               spearman_fisher_vs_logit = suppressWarnings(
                 cor(-log10(p_fisher), -log10(p_logit), method = "spearman")),
               pearson_beta_vs_log_or = suppressWarnings(
                 cor(beta, log(or_haldane), use = "complete.obs")),
               n_sig_logit_only  = sum(q_logit < FDR_ALPHA & !(q_perm_z < FDR_ALPHA), na.rm = TRUE),
               n_sig_perm_z_only = sum(q_perm_z < FDR_ALPHA & !(q_logit < FDR_ALPHA), na.rm = TRUE),
               n_sig_both        = sum(q_logit < FDR_ALPHA & q_perm_z < FDR_ALPHA, na.rm = TRUE)),
             by = .(gene_set, annotation_source)]

# =============================================================================
# 6. HOST-HOMOLOGY GRADIENT
#    Question: is evolutionary constraint associated with HAVING A DETECTABLE
#    SEQUENCE-SIMILAR human protein? This is sequence similarity only. It says
#    nothing about structural cross-reactivity, binding-site conservation or
#    actual off-target liability, none of which this analysis can establish.
# =============================================================================
msg("== 6. host-homology gradient ==")
hh <- copy(af_e)   # full constraint-evaluable universe; no annotation requirement
hh_rows <- list()

for (hv in c("sensitive", "strict")) {
  hh[, hhom := if (hv == "sensitive") has_human_homolog_sens else has_human_homolog_strict]

  # (i) constraint distribution by homolog status
  w <- wilcox.test(piNpiS_shrunk ~ hhom, data = hh)
  hh_rows[[length(hh_rows) + 1L]] <- data.table(
    test = "piNpiS by human-homolog status (Wilcoxon rank-sum)",
    homology_definition = hv, n = nrow(hh),
    statistic_name = "median piNpiS (homolog) - median piNpiS (no homolog)",
    estimate = median(hh$piNpiS_shrunk[hh$hhom]) - median(hh$piNpiS_shrunk[!hh$hhom]),
    ci_lo = NA_real_, ci_hi = NA_real_, p_raw = w$p.value, p_matched_null = NA_real_,
    note = "negative value = genes WITH a human homolog are more constrained")

  # (ii) unadjusted 2x2
  ft <- fisher.test(table(factor(hh$constr10, c(FALSE, TRUE)), factor(hh$hhom, c(FALSE, TRUE))))
  hh_rows[[length(hh_rows) + 1L]] <- data.table(
    test = "constrained(q10) x has-human-homolog, Fisher exact (unadjusted)",
    homology_definition = hv, n = nrow(hh), statistic_name = "odds ratio",
    estimate = unname(ft$estimate), ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2],
    p_raw = ft$p.value, p_matched_null = NA_real_, note = "")

  # (iii) logistic regression with covariates
  fit <- glm(hhom ~ constr10 + log10(cds_bp) + n_ortho, data = hh, family = binomial())
  s <- summary(fit)$coefficients
  hh_rows[[length(hh_rows) + 1L]] <- data.table(
    test = "logistic: has-human-homolog ~ constrained(q10) + log10(CDS bp) + n detectable Aspergillus orthologues",
    homology_definition = hv, n = nrow(hh), statistic_name = "beta(constrained) [log-odds]",
    estimate = s["constr10TRUE", 1], ci_lo = s["constr10TRUE", 1] - 1.96 * s["constr10TRUE", 2],
    ci_hi = s["constr10TRUE", 1] + 1.96 * s["constr10TRUE", 2],
    p_raw = s["constr10TRUE", 4], p_matched_null = NA_real_, note = "")

  # (iv) continuous gradient: constraint rank vs homolog indicator
  sp <- suppressWarnings(cor.test(rank(hh$piNpiS_shrunk), as.numeric(hh$hhom), method = "spearman"))
  hh_rows[[length(hh_rows) + 1L]] <- data.table(
    test = "Spearman: piNpiS rank vs has-human-homolog (continuous gradient)",
    homology_definition = hv, n = nrow(hh), statistic_name = "rho",
    estimate = unname(sp$estimate), ci_lo = NA_real_, ci_hi = NA_real_,
    p_raw = sp$p.value, p_matched_null = NA_real_,
    note = "negative rho = lower piNpiS (more constrained) associates with having a human homolog")

  # (v) MATCHED NULL: permute the constrained label within strata of
  #     (log10 CDS-length decile) x (conservation bin = n detectable Aspergillus
  #     orthologues, 0-6 / 7-8 / 9), statistic = fraction of the constrained set
  #     that has a detectable human homolog
  lg <- log10(hh$cds_bp)
  qs <- unique(quantile(lg, seq(0, 1, 0.1)))
  lbin <- cut(lg, qs, include.lowest = TRUE, labels = FALSE)
  cbin <- ifelse(hh$n_ortho <= 6, 1L, ifelse(hh$n_ortho <= 8, 2L, 3L))
  strat <- paste0("L", lbin, "_C", cbin)
  obs_frac <- mean(hh$hhom[hh$constr10])
  draw_hh <- make_sampler(strat, hh$constr10)
  nullf <- replicate(NPERM, mean(hh$hhom[draw_hh()]))
  p_hi <- (1 + sum(nullf >= obs_frac)) / (NPERM + 1)
  p_lo <- (1 + sum(nullf <= obs_frac)) / (NPERM + 1)
  p_z_hh <- 2 * pnorm(-abs((obs_frac - mean(nullf)) / sd(nullf)))
  hh_rows[[length(hh_rows) + 1L]] <- data.table(
    test = sprintf("matched-null permutation (%d perms; strata = log10-length decile x conservation bin): fraction of constrained genes with a human homolog",
                   NPERM),
    homology_definition = hv, n = nrow(hh), statistic_name = "observed fraction (null mean)",
    estimate = obs_frac, ci_lo = quantile(nullf, 0.025), ci_hi = quantile(nullf, 0.975),
    p_raw = NA_real_, p_matched_null = min(1, 2 * min(p_hi, p_lo)),
    note = sprintf("null mean %.4f, null sd %.4f, z = %.1f, z-based p = %.3g; CI columns are the 2.5-97.5%% matched-null interval; empirical p floor = %.2g",
                   mean(nullf), sd(nullf), (obs_frac - mean(nullf)) / sd(nullf),
                   p_z_hh, 2 / (NPERM + 1)))
  if (hv == "sensitive") { nullf_sens <- nullf; obs_frac_sens <- obs_frac }
}
hh_tab <- rbindlist(hh_rows)
hh_tab[, p_adj_BH := p.adjust(ifelse(is.na(p_raw), p_matched_null, p_raw), "BH")]
print(hh_tab[, .(test = substr(test, 1, 58), homology_definition,
                 estimate = round(estimate, 4), p = signif(fifelse(is.na(p_raw), p_matched_null, p_raw), 3),
                 q = signif(p_adj_BH, 3))])

# gradient table: % with a human homolog across constraint deciles
hh[, constraint_decile := cut(piNpiS_shrunk,
       breaks = quantile(piNpiS_shrunk, seq(0, 1, 0.1), na.rm = TRUE),
       include.lowest = TRUE, labels = 1:10)]
grad <- hh[, .(n = .N, median_piNpiS = median(piNpiS_shrunk),
               median_cds_bp = median(cds_bp), mean_n_ortho = mean(n_ortho),
               frac_human_homolog_sensitive = mean(has_human_homolog_sens),
               frac_human_homolog_strict = mean(has_human_homolog_strict),
               median_human_pident = median(human_pident, na.rm = TRUE)),
           by = constraint_decile][order(constraint_decile)]
grad[, decile_note := "1 = most constrained"]
print(grad)

# =============================================================================
# 7. OUTPUTS
# =============================================================================
msg("== 7. writing outputs ==")
setcolorder(res, c("gene_set", "annotation_source", "term", "term_name", "direction",
                   "n_universe", "n_universe_in_term", "n_set", "n_set_in_term",
                   "expected_in_term_hypergeom", "fold_enrichment",
                   "or_fisher", "or_lo", "or_hi", "or_haldane", "p_fisher", "q_fisher",
                   "beta", "beta_lo", "beta_hi", "p_logit", "q_logit", "logit_note",
                   "perm_null_mean", "perm_null_sd", "perm_fold", "perm_z",
                   "p_perm_enrich", "p_perm_deplete", "p_perm_two_sided", "q_perm",
                   "p_perm_z", "q_perm_z", "n_perm", "perm_empirical_floor",
                   "perm_at_empirical_floor", "q_fisher_global", "q_logit_global",
                   "q_perm_global", "q_perm_z_global",
                   "significant_any", "significant_all_three"))
setorder(res, gene_set, annotation_source, p_fisher)
fwrite(res, file.path(TAB, paste0("rev_b3_enrichment_all_tests", OUT_SUFFIX, ".csv")))

fwrite(univ_tab,      file.path(TAB, paste0("rev_b3_universe_definition", OUT_SUFFIX, ".csv")))
fwrite(gene_sets_tab, file.path(TAB, paste0("rev_b3_gene_sets", OUT_SUFFIX, ".csv")))
fwrite(ov_tab,        file.path(TAB, paste0("rev_b3_crossgenus_constraint_overlap", OUT_SUFFIX, ".csv")))
fwrite(sig_counts,    file.path(TAB, paste0("rev_b3_test_family_summary", OUT_SUFFIX, ".csv")))
fwrite(agree,         file.path(TAB, paste0("rev_b3_method_agreement", OUT_SUFFIX, ".csv")))
fwrite(match_tab,     file.path(TAB, paste0("rev_b3_permutation_matching_quality", OUT_SUFFIX, ".csv")))
fwrite(hh_tab,        file.path(TAB, paste0("rev_b3_host_homology_gradient", OUT_SUFFIX, ".csv")))
fwrite(grad,          file.path(TAB, paste0("rev_b3_host_homology_by_constraint_decile", OUT_SUFFIX, ".csv")))

sig <- res[significant_any == TRUE][order(gene_set, annotation_source, p_fisher)]
if (!nrow(sig)) {
  sig <- data.table(gene_set = "ALL", annotation_source = "ALL", term = NA_character_,
                    term_name = "NONE SURVIVED BH FDR 0.05 BY ANY OF THE THREE METHODS",
                    note = sprintf("%d tests across %d gene-set x annotation-source families",
                                   nrow(res), uniqueN(res[, .(gene_set, annotation_source)])))
}
fwrite(sig, file.path(TAB, paste0("rev_b3_significant_after_FDR", OUT_SUFFIX, ".csv")))

neg_sum <- res[grepl("^NEGCTRL", gene_set),
               .(n_terms = .N, n_sig_fisher_BH = sum(q_fisher < FDR_ALPHA, na.rm = TRUE),
                 n_sig_logit_BH = sum(q_logit < FDR_ALPHA, na.rm = TRUE),
                 n_sig_perm_empirical_BH = sum(q_perm < FDR_ALPHA, na.rm = TRUE),
                 n_sig_perm_z_BH = sum(q_perm_z < FDR_ALPHA, na.rm = TRUE),
                 frac_p_fisher_below_0.05 = mean(p_fisher < 0.05),
                 median_p_fisher = median(p_fisher),
                 ks_p_vs_uniform = suppressWarnings(ks.test(p_fisher, "punif")$p.value)),
               by = .(gene_set, annotation_source)]
fwrite(neg_sum,   file.path(TAB, paste0("rev_b3_negative_control_calibration", OUT_SUFFIX, ".csv")))
fwrite(calib_tab, file.path(TAB, paste0("rev_b3_negative_control_replicates", OUT_SUFFIX, ".csv")))
fwrite(calib_sum, file.path(TAB, paste0("rev_b3_negative_control_replicate_summary", OUT_SUFFIX, ".csv")))
print(neg_sum)

# ---------------------------- FIGURE ----------------------------------------
pdf(file.path(FIG, "rev_b3_pathway_enrichment.pdf"), width = 11, height = 13)
par(mfrow = c(3, 2), mar = c(4.3, 4.3, 3.0, 1.0), cex.axis = 0.9, cex.lab = 1.0)
pal <- c(AF_constrained_q10 = "#1f6fb4", CT_constrained_q10 = "#c1442e",
         BOTH_constrained_q25 = "#4b8b3b", NEGCTRL_AF_matched_random = "grey35",
         NEGCTRL_AF_uniform_random = "grey65")
show_sets <- names(pal)

# (A) p-value QQ against uniform -- calibration
plot(NA, xlim = c(0, 3.2), ylim = c(0, 3.2), xlab = expression(-log[10]~"expected p (uniform)"),
     ylab = expression(-log[10]~"observed p"), main = "A. Fisher p-values vs uniform (KEGG + GO)")
abline(0, 1, col = "grey70", lty = 2)
for (s in show_sets) {
  p <- sort(res[gene_set == s & annotation_source %in% c("KEGG_PATHWAY", "GO")]$p_fisher)
  if (!length(p)) next
  points(-log10(ppoints(length(p))), -log10(p), pch = 20, cex = 0.5, col = pal[[s]])
}
legend("topleft", legend = sub("_", " ", show_sets), col = pal, pch = 20, bty = "n", cex = 0.8)

# (B) volcano, constrained sets
v <- res[gene_set %in% c("AF_constrained_q10", "CT_constrained_q10", "BOTH_constrained_q25") &
           annotation_source %in% c("KEGG_PATHWAY", "GO")]
v[, l2or := log2(pmax(pmin(or_haldane, 64), 1/64))]
plot(v$l2or, -log10(v$p_fisher), pch = 20, cex = 0.55, col = pal[v$gene_set],
     xlab = expression(log[2]~"odds ratio (Haldane-Anscombe)"),
     ylab = expression(-log[10]~"p (Fisher)"),
     main = "B. Effect size vs significance")
abline(v = 0, col = "grey70", lty = 2)
thr <- v[q_fisher < FDR_ALPHA, if (.N) max(p_fisher) else NA_real_, by = gene_set]$V1
if (any(!is.na(thr))) abline(h = -log10(max(thr, na.rm = TRUE)), col = "firebrick", lty = 3)
lab <- v[order(p_fisher)][1:min(8, .N)]
text(lab$l2or, -log10(lab$p_fisher), substr(lab$term_name, 1, 26), cex = 0.52, pos = 4, xpd = NA)
legend("topleft", legend = c("BH q<0.05 boundary (if any)"), lty = 3, col = "firebrick",
       bty = "n", cex = 0.75)

# (C) logistic vs permutation agreement
a <- res[!is.na(p_logit) & !is.na(p_perm_two_sided) &
           annotation_source %in% c("KEGG_PATHWAY", "GO") & gene_set %in% show_sets]
rho_c <- if (nrow(a) > 3) suppressWarnings(cor(-log10(a$p_logit),
            -log10(a$p_perm_two_sided), method = "spearman")) else NA_real_
plot(-log10(a$p_logit), -log10(a$p_perm_two_sided), pch = 20, cex = 0.55,
     col = pal[a$gene_set], xlab = expression(-log[10]~"p (logistic, covariate-adjusted)"),
     ylab = expression(-log[10]~"p (matched permutation, empirical)"),
     main = sprintf("C. Method agreement (Spearman rho = %.2f)", rho_c))
abline(0, 1, col = "grey70", lty = 2)
abline(h = -log10(2 / (NPERM + 1)), col = "firebrick", lty = 3)
legend("topleft", legend = sprintf("permutation resolution floor (2/%d)", NPERM + 1),
       lty = 3, col = "firebrick", bty = "n", cex = 0.72)

# (D) constraint distribution by human-homolog status
b1 <- hh$piNpiS_shrunk[!hh$has_human_homolog_sens]
b2 <- hh$piNpiS_shrunk[hh$has_human_homolog_sens]
boxplot(list(`no human\nhomolog` = b1, `human homolog\n(sequence similarity)` = b2),
        outline = FALSE, col = c("grey85", "#1f6fb4"), ylab = "piN/piS (shrunk)",
        main = "D. Constraint by human sequence similarity")
mtext(sprintf("medians %.3f vs %.3f", median(b1), median(b2)), side = 3, cex = 0.72, line = 0.1)

# (E) matched-null for the host-homology statistic
hist(nullf_sens, breaks = 40, col = "grey85", border = "white",
     xlim = range(c(nullf_sens, obs_frac_sens)) + c(-0.012, 0.012),
     xlab = "fraction of the constrained set with a human homolog",
     main = "E. Host-homology vs length/conservation-matched null")
abline(v = obs_frac_sens, col = "firebrick", lwd = 2)
legend("topright", legend = c(sprintf("observed = %.3f", obs_frac_sens),
                              sprintf("matched null = %.3f", mean(nullf_sens))),
       col = c("firebrick", "grey60"), lwd = c(2, 6), bty = "n", cex = 0.78)

# (F) gradient across constraint deciles
par(mar = c(6.0, 4.3, 3.0, 1.0))
plot(as.integer(grad$constraint_decile), grad$frac_human_homolog_sensitive, type = "b",
     pch = 19, col = "#1f6fb4", ylim = c(0, 1), xlab = "piN/piS decile (1 = most constrained)",
     ylab = "fraction with a detectable human homolog",
     main = "F. Host-similarity gradient across constraint deciles")
lines(as.integer(grad$constraint_decile), grad$frac_human_homolog_strict, type = "b",
      pch = 17, col = "#c1442e")
legend("bottomleft", legend = c("sensitive BLAST set", "strict BLAST set"),
       col = c("#1f6fb4", "#c1442e"), pch = c(19, 17), bty = "n", cex = 0.8)
mtext("sequence similarity only; NOT evidence of structural", side = 1, line = 3.9, cex = 0.62)
mtext("cross-reactivity or off-target risk", side = 1, line = 4.7, cex = 0.62)
dev.off()

msg("== done ==")
msg("tests written: %d rows over %d gene-set x annotation-source families",
    nrow(res), uniqueN(res[, .(gene_set, annotation_source)]))
msg("terms surviving BH %.2f by any method: %d",
    FDR_ALPHA, sum(res$significant_any, na.rm = TRUE))
