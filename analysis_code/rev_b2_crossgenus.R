#!/usr/bin/env Rscript
# rev_b2: does the published cross-genus claim ("constraint does not transfer between
# genera, Spearman rho ~ 0.27") survive matching the estimator and the sampling
# between A. fumigatus and C. tropicalis?
#
# The published comparison (analysis_code/crossgenus_intersection.R) is MISMATCHED on
# three axes at once:
#   estimator : Af uses the pseudocount-shrunk count ratio, Ct the raw pi ratio
#   sampling  : Af uses all 315 genomes, Ct uses 151 lineage representatives
#   depth     : median segregating sites per gene ~14 (Af) vs ~108 (Ct)
# Here each axis is matched in turn and in combination. Within-genus split-half
# reliability gives the ceiling that estimation noise alone imposes on rho, so the
# residual (disattenuated) rho is the biological transferability.
#
# Outputs: output/table/rev_b2_crossgenus_matched.csv
#          output/table/rev_b2_crossgenus_depthmatch.csv
#          output/table/rev_b2_reliability.csv
suppressPackageStartupMessages({library(data.table); library(parallel)})
setDTthreads(8)
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
TD <- "analysis/target_discovery"; CAN <- "analysis/candida"; OT <- "output/table"
SD <- Sys.getenv("REV_B2_SCRATCH",
  "/tmp/claude-1000/-mnt-HuangLab-26T-HuangLab-Data-Tar-Nine-Ta/459df96a-3261-4f14-8f94-382b2767e924/scratchpad")
NBOOT <- 2000; NCORE <- 24
set.seed(20261012)

af_stat <- fread(file.path(TD, "rev_b2_gene_static.tsv"))
af_cnt  <- fread(cmd = paste("zcat", file.path(TD, "rev_b2_setting_counts.tsv.gz")))
ct_stat <- fread(file.path(CAN, "constraint/rev_b2_ctrop_gene_static.tsv"))
ct_cnt  <- fread(cmd = paste("zcat", file.path(CAN, "constraint/rev_b2_ctrop_setting_counts.tsv.gz")))
om      <- fread("analysis/crossgenus/ortho_map.tsv")
pub_af  <- fread(file.path(TD, "constraint_per_gene.tsv"))
pub_ct  <- fread(file.path(CAN, "constraint/ctrop_constraint_per_gene.tsv"))

## ------------------------------------------------- estimators (identical formulae)
est <- function(d, kind) {
  if (kind == "raw")        fifelse(d$pi_syn > 0, d$pi_nonsyn / d$pi_syn, NA_real_)
  else if (kind == "genomemean") {
    gN <- mean(d$pi_nonsyn, na.rm = TRUE); gS <- mean(d$pi_syn, na.rm = TRUE)
    (d$pi_nonsyn + gN) / (d$pi_syn + gS)
  } else {
    pc <- as.numeric(sub("^pc", "", kind))
    ((d$nonsyn_seg + pc) / d$nonsyn_sites) / ((d$syn_seg + pc) / d$syn_sites)
  }
}

af_tab <- function(sname, kind, callab_frac = 250 / 315) {
  x <- af_cnt[setting == sname]
  n <- x$n_genomes[1]
  d <- merge(af_stat[, .(protein_id, syn_sites, nonsyn_sites)], x, by = "protein_id")
  d[, cstat := est(d, kind)]
  d[mean_cov_genomes >= callab_frac * n & (syn_seg + nonsyn_seg) >= 3 & is.finite(cstat),
    .(af_protein = protein_id, af_pi = cstat, af_seg = syn_seg + nonsyn_seg)]
}
ct_tab <- function(sname, kind, min_call = 300) {
  x <- ct_cnt[setting == sname]
  d <- merge(ct_stat[, .(protein_id, syn_sites, nonsyn_sites, n_iso_callable)], x, by = "protein_id")
  d[, cstat := est(d, kind)]
  d[n_iso_callable >= min_call & (syn_seg + nonsyn_seg) >= 3 & is.finite(cstat),
    .(ct_protein = protein_id, ct_pi = cstat, ct_seg = syn_seg + nonsyn_seg)]
}

pairs_base <- om[ct_locus != "", .(af_protein, af_locus, ct_protein, ct_locus)]

join_pairs <- function(a, c) {
  m <- merge(pairs_base, a, by = "af_protein")
  merge(m, c, by = "ct_protein")
}

boot_rho <- function(x, y, nboot = NBOOT) {
  n <- length(x)
  # The resampling indices are drawn here, in the parent process, from the one stream that
  # set.seed() above controls. Drawing them inside mclapply() instead is not reproducible:
  # with the default mc.set.seed = TRUE and Mersenne-Twister the workers are seeded from time
  # and pid, so the top-level seed never reaches them and the interval moves by a few parts
  # per thousand on every run. Drawing them in the parent also makes the interval independent
  # of how many cores the run is given. Only the correlation is parallelised.
  idx <- matrix(sample.int(n, n * nboot, replace = TRUE), nrow = n, ncol = nboot)
  r <- unlist(mclapply(seq_len(nboot), function(i) {
    j <- idx[, i]
    suppressWarnings(cor(x[j], y[j], method = "spearman"))
  }, mc.cores = NCORE))
  quantile(r, c(0.025, 0.975), na.rm = TRUE)
}

one <- function(label, af_set, af_kind, ct_set, ct_kind, note = "") {
  m <- join_pairs(af_tab(af_set, af_kind), ct_tab(ct_set, ct_kind))
  rho <- suppressWarnings(cor(m$af_pi, m$ct_pi, method = "spearman"))
  tau <- suppressWarnings(cor(m$af_pi, m$ct_pi, method = "kendall"))
  ci <- boot_rho(m$af_pi, m$ct_pi)
  p <- suppressWarnings(cor.test(m$af_pi, m$ct_pi, method = "spearman",
                                 exact = FALSE))$p.value
  data.table(comparison = label, af_setting = af_set, af_estimator = af_kind,
             ct_setting = ct_set, ct_estimator = ct_kind, note = note,
             n_pairs = nrow(m), spearman_rho = rho, rho_lo95 = ci[[1]], rho_hi95 = ci[[2]],
             kendall_tau = tau, p_value = p,
             af_median_seg = median(m$af_seg), ct_median_seg = median(m$ct_seg))
}

## ------------------------------------------------- 1. reproduce the published value
pa <- merge(af_stat[, .(protein_id, syn_sites, nonsyn_sites)],
            pub_af[, .(protein_id, syn_seg, nonsyn_seg, pi_syn, pi_nonsyn, mean_cov_genomes)],
            by = "protein_id")
pa[, cstat := est(pa, "pc0.5")]
pa <- pa[mean_cov_genomes >= 250 & (syn_seg + nonsyn_seg) >= 3 & is.finite(cstat),
         .(af_protein = protein_id, af_pi = cstat, af_seg = syn_seg + nonsyn_seg)]
pc_ <- pub_ct[n_iso_callable >= 300 & is.finite(piNpiS_dedup) & (syn_seg + nonsyn_seg) >= 3,
              .(ct_protein = protein_id, ct_pi = piNpiS_dedup, ct_seg = syn_seg + nonsyn_seg)]
mp <- join_pairs(pa, pc_)
rho_pub <- suppressWarnings(cor(mp$af_pi, mp$ct_pi, method = "spearman"))
ci_pub <- boot_rho(mp$af_pi, mp$ct_pi)
cat(sprintf("PUBLISHED comparison reproduced: n = %d pairs, Spearman rho = %.3f [%.3f, %.3f]\n",
            nrow(mp), rho_pub, ci_pub[[1]], ci_pub[[2]]))
cat(sprintf("  median segregating sites/gene: Af %.0f vs Ct %.0f\n",
            median(mp$af_seg), median(mp$ct_seg)))

R <- rbindlist(list(
  data.table(comparison = "PUBLISHED (mismatched estimator + sampling)",
             af_setting = "constraint_per_gene.tsv (all 315)", af_estimator = "pc0.5 shrunk",
             ct_setting = "ctrop piNpiS_dedup (151 reps)", ct_estimator = "raw",
             note = "as in crossgenus_intersection.R", n_pairs = nrow(mp),
             spearman_rho = rho_pub, rho_lo95 = ci_pub[[1]], rho_hi95 = ci_pub[[2]],
             kendall_tau = suppressWarnings(cor(mp$af_pi, mp$ct_pi, method = "kendall")),
             p_value = suppressWarnings(cor.test(mp$af_pi, mp$ct_pi, method = "spearman",
                                                 exact = FALSE))$p.value,
             af_median_seg = median(mp$af_seg), ct_median_seg = median(mp$ct_seg)),
  ## 2. match the ESTIMATOR only (sampling left as published)
  one("estimator matched: both raw", "all315_mac1", "raw", "ct_dedup_mac1", "raw",
      "Af sampling all 315, Ct lineage reps (as published)"),
  one("estimator matched: both pc0.5 shrunk", "all315_mac1", "pc0.5", "ct_dedup_mac1", "pc0.5",
      "Af sampling all 315, Ct lineage reps (as published)"),
  one("estimator matched: both genome-mean shrunk", "all315_mac1", "genomemean",
      "ct_dedup_mac1", "genomemean", "Af sampling all 315, Ct lineage reps"),
  ## 3. match the SAMPLING only (estimator left as published)
  one("sampling matched: both all genomes/isolates", "all315_mac1", "pc0.5", "ct_all_mac1", "raw",
      "estimator still mismatched"),
  one("sampling matched: both lineage-deduplicated", "dedup188_mac1", "pc0.5",
      "ct_dedup_mac1", "raw", "estimator still mismatched"),
  ## 4. match BOTH
  one("MATCHED: both raw, both all", "all315_mac1", "raw", "ct_all_mac1", "raw", "fully matched"),
  one("MATCHED: both raw, both dedup", "dedup188_mac1", "raw", "ct_dedup_mac1", "raw", "fully matched"),
  one("MATCHED: both pc0.5, both all", "all315_mac1", "pc0.5", "ct_all_mac1", "pc0.5", "fully matched"),
  one("MATCHED: both pc0.5, both dedup", "dedup188_mac1", "pc0.5", "ct_dedup_mac1", "pc0.5", "fully matched"),
  one("MATCHED: both genome-mean, both dedup", "dedup188_mac1", "genomemean",
      "ct_dedup_mac1", "genomemean", "fully matched"),
  one("MATCHED: both pc0.5, both dedup, both mac>=2", "dedup188_mac2", "pc0.5",
      "ct_dedup_mac2", "pc0.5", "fully matched + singletons dropped in both")
))
print(R[, .(comparison, n = n_pairs, rho = round(spearman_rho, 3),
            CI = sprintf("[%.3f, %.3f]", rho_lo95, rho_hi95), tau = round(kendall_tau, 3),
            af_seg = af_median_seg, ct_seg = ct_median_seg)])
fwrite(R, file.path(OT, "rev_b2_crossgenus_matched.csv"))

## ------------------------------------------------- 5. matched SAMPLING DEPTH
# Af has ~14 segregating sites per gene, Ct ~108. Shrink the Ct isolate set until the
# two genera carry comparable information per gene, then re-measure rho.
ds <- grep("^ct_dedupsub", unique(ct_cnt$setting), value = TRUE)
D <- rbindlist(mclapply(ds, function(s) {
  nsub <- as.integer(sub("^ct_dedupsub(\\d+)_r\\d+$", "\\1", s))
  out <- rbindlist(lapply(c("raw", "pc0.5"), function(kind) {
    m <- join_pairs(af_tab(if (kind == "raw") "dedup188_mac1" else "dedup188_mac1", kind),
                    ct_tab(s, kind))
    if (nrow(m) < 50) return(NULL)
    data.table(ct_setting = s, ct_n_iso = nsub, estimator = kind, n_pairs = nrow(m),
               spearman_rho = suppressWarnings(cor(m$af_pi, m$ct_pi, method = "spearman")),
               af_median_seg = median(m$af_seg), ct_median_seg = median(m$ct_seg))
  }))
  out
}, mc.cores = NCORE))
Dsum <- D[, .(n_draws = .N, n_pairs_median = as.numeric(median(as.numeric(n_pairs))),
              ct_median_seg = as.numeric(median(as.numeric(ct_median_seg))),
              af_median_seg = as.numeric(median(as.numeric(af_median_seg))),
              rho_mean = mean(spearman_rho), rho_sd = sd(spearman_rho),
              rho_min = min(spearman_rho), rho_max = max(spearman_rho)),
          by = .(ct_n_iso, estimator)][order(estimator, ct_n_iso)]
cat("\n=== matched sampling depth (both genera lineage-deduplicated, matched estimator) ===\n")
print(Dsum)
fwrite(D, file.path(OT, "rev_b2_crossgenus_depthmatch.csv"))
fwrite(Dsum, file.path(OT, "rev_b2_crossgenus_depthmatch_summary.csv"))

## ------------------------------------------------- 6. within-genus reliability
# Split-half Spearman between two disjoint halves of each genus' lineage representatives.
# Spearman-Brown corrects a half-sample correlation up to the full-sample reliability.
afsp <- fread(cmd = paste("zcat", file.path(TD, "rev_b2_splithalf_counts.tsv.gz")))
ctsp <- fread(cmd = paste("zcat", file.path(CAN, "constraint/rev_b2_ctrop_splithalf_counts.tsv.gz")))

rel_af <- function(kind) {
  dr <- unique(sub("[AB]$", "", afsp$setting))
  sapply(dr, function(r) {
    a <- merge(af_stat[, .(protein_id, syn_sites, nonsyn_sites)], afsp[setting == paste0(r, "A")], by = "protein_id")
    b <- merge(af_stat[, .(protein_id, syn_sites, nonsyn_sites)], afsp[setting == paste0(r, "B")], by = "protein_id")
    a[, cstat := est(a, kind)]; b[, cstat := est(b, kind)]
    n <- a$n_genomes[1]
    ka <- a[mean_cov_genomes >= (250/315) * n & (syn_seg + nonsyn_seg) >= 3 & is.finite(cstat)]
    kb <- b[mean_cov_genomes >= (250/315) * n & (syn_seg + nonsyn_seg) >= 3 & is.finite(cstat)]
    mm <- merge(ka[, .(protein_id, x = cstat)], kb[, .(protein_id, y = cstat)], by = "protein_id")
    suppressWarnings(cor(mm$x, mm$y, method = "spearman"))
  })
}
rel_ct <- function(kind) {
  dr <- unique(sub("[AB]$", "", ctsp$setting))
  sapply(dr, function(r) {
    a <- merge(ct_stat[, .(protein_id, syn_sites, nonsyn_sites, n_iso_callable)], ctsp[setting == paste0(r, "A")], by = "protein_id")
    b <- merge(ct_stat[, .(protein_id, syn_sites, nonsyn_sites, n_iso_callable)], ctsp[setting == paste0(r, "B")], by = "protein_id")
    a[, cstat := est(a, kind)]; b[, cstat := est(b, kind)]
    ka <- a[n_iso_callable >= 300 & (syn_seg + nonsyn_seg) >= 3 & is.finite(cstat)]
    kb <- b[n_iso_callable >= 300 & (syn_seg + nonsyn_seg) >= 3 & is.finite(cstat)]
    mm <- merge(ka[, .(protein_id, x = cstat)], kb[, .(protein_id, y = cstat)], by = "protein_id")
    suppressWarnings(cor(mm$x, mm$y, method = "spearman"))
  })
}
sb <- function(r) 2 * r / (1 + r)       # Spearman-Brown half -> full
REL <- rbindlist(lapply(c("raw", "pc0.5"), function(k) {
  a <- rel_af(k); c_ <- rel_ct(k)
  data.table(estimator = k,
             af_splithalf_rho = mean(a), af_splithalf_sd = sd(a),
             af_reliability_SB = sb(mean(a)),
             ct_splithalf_rho = mean(c_), ct_splithalf_sd = sd(c_),
             ct_reliability_SB = sb(mean(c_)), n_splits = length(a))
}))
cat("\n=== within-genus split-half reliability (disjoint halves of lineage representatives) ===\n")
print(REL)

## disattenuated cross-genus rho for the fully matched comparisons
dis <- R[grepl("^MATCHED", comparison)]
dis[, estimator := fifelse(grepl("raw", comparison), "raw",
                    fifelse(grepl("pc0.5", comparison), "pc0.5", NA_character_))]
dis <- merge(dis[, .(comparison, estimator, n_pairs, spearman_rho, rho_lo95, rho_hi95)],
             REL[, .(estimator, af_reliability_SB, ct_reliability_SB)], by = "estimator", all.x = TRUE)
dis[, rho_disattenuated := spearman_rho / sqrt(af_reliability_SB * ct_reliability_SB)]
dis[, rho_disattenuated_lo := rho_lo95 / sqrt(af_reliability_SB * ct_reliability_SB)]
dis[, rho_disattenuated_hi := rho_hi95 / sqrt(af_reliability_SB * ct_reliability_SB)]
cat("\n=== noise-corrected cross-genus rho ===\n")
print(dis[, .(comparison, rho = round(spearman_rho, 3),
              rel_Af = round(af_reliability_SB, 3), rel_Ct = round(ct_reliability_SB, 3),
              rho_corrected = round(rho_disattenuated, 3),
              corrected_CI = sprintf("[%.3f, %.3f]", rho_disattenuated_lo, rho_disattenuated_hi))])
fwrite(merge(REL, unique(dis[, .(estimator)]), by = "estimator", all = TRUE),
       file.path(OT, "rev_b2_reliability.csv"))
fwrite(dis, file.path(OT, "rev_b2_crossgenus_disattenuated.csv"))

saveRDS(list(R = R, D = D, Dsum = Dsum, REL = REL, dis = dis, rho_pub = rho_pub,
             ci_pub = ci_pub, n_pub = nrow(mp)),
        file.path(SD, "rev_b2_crossgenus_results.rds"))
cat("\ndone\n")
