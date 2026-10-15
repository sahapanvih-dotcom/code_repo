#!/usr/bin/env Rscript
# rev_b2: robustness / sensitivity of the A. fumigatus evolutionary-constraint ranking.
#
# Re-runs the published funnel (analysis_code/td_funnel.R) over a grid of analytical
# settings, using per-gene variant counts recomputed by
# analysis_code/rev_b2_asp_recompute.py (which reuses the project's existing
# td_constraint.py parser; everything downstream of the counts is done here in R).
#
# Outputs: output/table/rev_b2_settings_metrics.csv
#          output/table/rev_b2_resample_metrics.csv
#          output/table/rev_b2_tier1_retention.csv
#          output/table/rev_b2_candidate_percentiles.csv
#          output/table/rev_b2_funnel_counts.csv
suppressPackageStartupMessages({
  library(data.table); library(parallel); library(pcaPP)
})
setDTthreads(8)
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
TD <- "analysis/target_discovery"
OT <- "output/table"
NCORE <- 40

## ---------------------------------------------------------------- static data
stat <- fread(file.path(TD, "rev_b2_gene_static.tsv"))
cnt  <- fread(cmd = paste("zcat", file.path(TD, "rev_b2_setting_counts.tsv.gz")))

## selectivity / evidence layers - identical to td_funnel.R, setting-independent
hum <- fread(file.path(TD, "hom/af293_vs_human.tsv"),
             col.names = c("protein_id","human_hit","pident","len","evalue","bits","qcov","scov"))
hum <- hum[, .SD[which.min(evalue)], by = protein_id]
par_ <- fread(file.path(TD, "hom/has_paralog.txt"), header = FALSE)$V1
xs  <- fread(file.path(TD, "hom/crossspecies_besthits.tsv"), header = FALSE,
             col.names = c("protein_id","species"))
xsn <- xs[, .(n_species = uniqueN(species)), by = protein_id]
known <- c("AFUA_4G06890","AFUA_7G03740","AFUA_6G12400","AFUA_2G11010","AFUA_1G14870",
           "AFUA_2G03700","AFUA_1G14330","AFUA_2G02690","AFUA_5G04170","AFUA_5G09360",
           "AFUA_4G04820","AFUA_6G05300","AFUA_2G01260","AFUA_3G09850")

layers <- merge(stat, hum[, .(protein_id, human_evalue = evalue)], by = "protein_id", all.x = TRUE)
layers <- merge(layers, xsn, by = "protein_id", all.x = TRUE)
layers[is.na(n_species), n_species := 0L]
layers[, has_paralog := protein_id %in% par_]
layers[, has_human_homolog := !is.na(human_evalue)]
layers[, known_target := locus_tag %in% known]
layers[, cds_bp := nonsyn_sites + syn_sites]

## published reference ranking + shortlist
pub_raw  <- fread(file.path(TD, "constraint_per_gene.tsv"))
pub_cand <- fread("doc/Aspergillus/target_discovery_2026-08/tables/candidate_targets_ranked.csv")
tier1    <- fread(file.path(OT, "ms_table3_tier1.csv"))
N_GENOMES_PUB <- 315

## ------------------------------------------------------- the funnel, generalised
# One function implementing td_funnel.R, with every published threshold exposed.
# Baseline arguments reproduce td_funnel.R exactly.
run_funnel <- function(dc, n_gen,
                       callab_frac = 250 / 315,   # mean_cov_genomes >= 250 of 315
                       shrink      = "pc0.5",     # pc<x> | none | genomemean
                       cutoff      = "pub",       # pub|strict|loose|none|top10pct|top25pct
                       power_min   = 15,          # exp_tot >= 15
                       power_mode  = "exp_tot",   # exp_tot | cds_bp (length-matched)
                       dens_ref    = NA_real_,    # baseline variant density, for cds_bp mode
                       eval_rule   = "pub") {     # pub | seg5
  d <- copy(dc)
  d[, tot_seg := syn_seg + nonsyn_seg]
  d[, exp_nonsyn_frac := nonsyn_sites / (nonsyn_sites + syn_sites)]
  d[, p_nonsyn_depl := ifelse(tot_seg >= 3, pbinom(nonsyn_seg, tot_seg, exp_nonsyn_frac), NA_real_)]
  dens <- sum(d$tot_seg) / sum(d$cds_bp)
  d[, exp_tot := cds_bp * dens]
  d[, p_tot_depl := ppois(tot_seg, exp_tot)]

  ## constraint statistic (lower = more constrained)
  if (shrink == "none") {
    d[, cstat := fifelse(pi_syn > 0, pi_nonsyn / pi_syn, NA_real_)]
  } else if (shrink == "genomemean") {
    gN <- mean(d$pi_nonsyn, na.rm = TRUE); gS <- mean(d$pi_syn, na.rm = TRUE)
    d[, cstat := (pi_nonsyn + gN) / (pi_syn + gS)]
  } else {
    pc <- as.numeric(sub("^pc", "", shrink))
    d[, cstat := ((nonsyn_seg + pc) / nonsyn_sites) / ((syn_seg + pc) / syn_sites)]
  }

  d[, callable := mean_cov_genomes >= callab_frac * n_gen]
  d[, evaluable := if (eval_rule == "seg5") tot_seg >= 5 else (tot_seg >= 3 | exp_tot >= 10)]

  ref <- d[callable == TRUE & tot_seg >= 3]
  rd  <- ref$cstat
  # percentile against the callable+evaluable reference distribution, as in td_funnel.R
  d[, constraint_pctile := 100 * sapply(cstat, function(x) if (is.na(x)) NA_real_ else mean(rd > x, na.rm = TRUE))]

  f1 <- d[callable == TRUE]
  f2 <- f1[evaluable == TRUE]
  if (cutoff == "pub")           f3 <- f2[p_nonsyn_depl < 0.05 | p_tot_depl < 0.01]
  else if (cutoff == "strict")   f3 <- f2[p_nonsyn_depl < 0.01 | p_tot_depl < 0.001]
  else if (cutoff == "loose")    f3 <- f2[p_nonsyn_depl < 0.10 | p_tot_depl < 0.05]
  else if (cutoff == "none")     f3 <- f2
  else if (cutoff == "top10pct") f3 <- f2[cstat <= quantile(f2$cstat, 0.10, na.rm = TRUE)]
  else if (cutoff == "top25pct") f3 <- f2[cstat <= quantile(f2$cstat, 0.25, na.rm = TRUE)]
  else stop("cutoff")
  f3 <- f3[!is.na(cstat)]
  f4 <- f3[has_human_homolog == FALSE]
  f5 <- f4[has_paralog == FALSE]
  f6 <- f5[n_species >= 8]
  f7 <- f6[known_target == FALSE]
  f8 <- f7[protein_len >= 100]
  if (power_mode == "cds_bp") {
    # exp_tot = cds_bp * dens, so an absolute exp_tot gate moves whenever a
    # variant filter changes the genome-wide density. This mode holds the gate at
    # the SAME SET OF GENES as the baseline (gate on gene length at baseline
    # density), isolating the effect of the variant filter from the gate shift.
    f9 <- f8[cds_bp >= power_min / dens_ref]
  } else {
    f9 <- f8[exp_tot >= power_min]
  }
  f9[, constraint_score := -log10(pmin(ifelse(is.na(p_nonsyn_depl), 1, p_nonsyn_depl),
                                       ifelse(is.na(p_tot_depl), 1, p_tot_depl)) + 1e-300)]
  setorder(f9, cstat, -constraint_score)
  f9[, rank := .I]
  list(all = d,
       dens = dens,
       funnel = c(all = nrow(d), callable = nrow(f1), evaluable = nrow(f2), constrained = nrow(f3),
                  no_human = nrow(f4), no_paralog = nrow(f5), conserved = nrow(f6),
                  not_known = nrow(f7), size_ok = nrow(f8), powered = nrow(f9)),
       ranked = f9)
}

jac <- function(a, b, k) {
  A <- head(a, k); Bb <- head(b, k)
  length(intersect(A, Bb)) / length(union(A, Bb))
}

## ---------------------------------------------------------------- baseline
mk <- function(sname) {
  x <- cnt[setting == sname]
  n <- x$n_genomes[1]
  d <- merge(layers, x[, .(protein_id, syn_seg, nonsyn_seg, stop_seg, pi_syn, pi_nonsyn,
                           mean_cov_genomes)], by = "protein_id")
  list(d = d, n = n)
}

base <- mk("all315_mac1")
pub <- run_funnel(base$d, base$n)
DENS_REF <- pub$dens
cat("=== baseline funnel reproduction (should match td_funnel.R) ===\n")
cat(sprintf("baseline genome-wide variant density = %.6g per CDS bp\n", DENS_REF))
print(pub$funnel)
cat("published candidate_targets_ranked.csv rows:", nrow(pub_cand), "\n")
cat("my baseline ranked rows:", nrow(pub$ranked), "\n")
cat("locus overlap with published shortlist:",
    length(intersect(pub$ranked$locus_tag, pub_cand$locus_tag)), "\n")
cat("top-20 Jaccard vs published order:",
    round(jac(pub$ranked$locus_tag, pub_cand$locus_tag, 20), 3), "\n\n")

## reference ranking used for all comparisons = the PUBLISHED file
pub_all <- run_funnel(merge(layers, pub_raw[, .(protein_id, syn_seg, nonsyn_seg, stop_seg,
                                                pi_syn, pi_nonsyn, mean_cov_genomes)],
                            by = "protein_id"), N_GENOMES_PUB)
ref_stat <- pub_all$all[, .(protein_id, locus_tag, ref_cstat = cstat,
                            ref_eval = callable & evaluable)]
ref_order <- pub_cand$locus_tag            # published ranked shortlist order
shortlist <- pub_cand$locus_tag
t1 <- tier1$locus_tag

## ---------------------------------------------------------------- metric engine
metrics <- function(res, label, factor_name, level) {
  d <- merge(res$all[, .(protein_id, locus_tag, cstat, callable, evaluable, constraint_pctile)],
             ref_stat, by = c("protein_id","locus_tag"))
  g <- d[ref_eval == TRUE & callable == TRUE & evaluable == TRUE & !is.na(cstat) & !is.na(ref_cstat)]
  sl <- g[locus_tag %in% shortlist]
  rk <- res$ranked$locus_tag
  data.table(
    setting = label, factor = factor_name, level = level,
    n_genes_compared = nrow(g),
    rho_genomewide = suppressWarnings(cor(g$cstat, g$ref_cstat, method = "spearman")),
    tau_genomewide = if (nrow(g) > 10) cor.fk(frank(g$cstat), frank(g$ref_cstat)) else NA_real_,
    n_shortlist_compared = nrow(sl),
    rho_shortlist = if (nrow(sl) > 10) suppressWarnings(cor(sl$cstat, sl$ref_cstat, method = "spearman")) else NA_real_,
    tau_shortlist = if (nrow(sl) > 10) cor.fk(frank(sl$cstat), frank(sl$ref_cstat)) else NA_real_,
    n_eligible = length(rk),
    # set-level Jaccard of the WHOLE eligible shortlist: separates "who is eligible"
    # from "how the eligible genes are ordered" (which is what top-k Jaccard mixes).
    jaccard_eligible_set = length(intersect(rk, ref_order)) / length(union(rk, ref_order)),
    frac_published_shortlist_retained = mean(ref_order %in% rk),
    jaccard_top20 = jac(rk, ref_order, 20),
    jaccard_top50 = jac(rk, ref_order, 50),
    jaccard_top100 = jac(rk, ref_order, 100),
    overlap_top20 = length(intersect(head(rk, 20), head(ref_order, 20))),
    overlap_top50 = length(intersect(head(rk, 50), head(ref_order, 50))),
    overlap_top100 = length(intersect(head(rk, 100), head(ref_order, 100))),
    tier1_eligible = sum(t1 %in% rk),
    tier1_in_top100 = sum(t1 %in% head(rk, 100)),
    tier1_in_top200 = sum(t1 %in% head(rk, 200)),
    top20_median_nonsyn_seg = if (length(rk) >= 20)
      as.numeric(median(res$ranked$nonsyn_seg[1:20])) else NA_real_)
}

## ---------------------------------------------------------------- setting grid
BASE_ARGS <- list(callab_frac = 250/315, shrink = "pc0.5", cutoff = "pub",
                  power_min = 15, power_mode = "exp_tot", evaluable = "pub")  # 'evaluable' -> eval_rule in run_funnel
grid <- list()
add <- function(label, data_setting, factor_name, level, ...) {
  a <- modifyList(BASE_ARGS, list(...))
  grid[[length(grid) + 1]] <<- c(list(label = label, data = data_setting,
                                      factor = factor_name, level = level), a)
}

add("baseline (published)", "all315_mac1", "baseline", "published")
## 1. derived-allele-count / frequency threshold
for (s in c("all315_mac2","all315_mac3","all315_mac5","all315_maf0.01","all315_maf0.05"))
  add(s, s, "allele count threshold", sub("all315_", "", s))
## 2. lineage de-duplication
add("dedup188_mac1", "dedup188_mac1", "population sampling", "lineage representatives (n=188)")
## 3. callability threshold
for (cf in c(200, 280, 300))
  add(sprintf("callab>=%d/315", cf), "all315_mac1", "callability threshold",
      sprintf("%d/315 (%.2f)", cf, cf/315), callab_frac = cf/315)
## 4. shrinkage
for (sh in c("none","pc0.25","pc1","pc2","genomemean"))
  add(paste0("shrink=", sh), "all315_mac1", "shrinkage estimator", sh, shrink = sh)
## 5. constrained-gene cutoff / top percentile
for (co in c("strict","loose","none","top10pct","top25pct"))
  add(paste0("cutoff=", co), "all315_mac1", "constrained cutoff", co, cutoff = co)
## 6. prespecified expected-variant power gate
for (pw in c(0, 10, 25, 40))
  add(sprintf("exp_tot>=%d", pw), "all315_mac1", "expected-variant gate",
      sprintf("exp_tot>=%d", pw), power_min = pw)
## 7. evaluable gate
add("evaluable=tot_seg>=5", "all315_mac1", "evaluable gate", "tot_seg>=5", evaluable = "seg5")
## 8. combinations
add("dedup + no shrinkage", "dedup188_mac1", "combination", "dedup + raw ratio", shrink = "none")
add("dedup + mac2", "dedup188_mac2", "combination", "dedup + mac2")
add("dedup + mac3", "dedup188_mac3", "combination", "dedup + mac3")
add("dedup + maf0.05", "dedup188_maf0.05", "combination", "dedup + maf0.05")
add("mac2 + pc1 + exp40", "all315_mac2", "combination", "mac2 + pc1 + exp_tot>=40",
    shrink = "pc1", power_min = 40)
add("strict cutoff + exp40 + callab300", "all315_mac1", "combination",
    "strict + exp_tot>=40 + callab 300/315", cutoff = "strict", power_min = 40,
    callab_frac = 300/315)
add("dedup + mac3 + raw + strict", "dedup188_mac3", "combination",
    "dedup + mac3 + raw + strict", shrink = "none", cutoff = "strict")
## 9. allele-count thresholds with the power gate held at the baseline gene set
for (s2 in c("all315_mac2","all315_mac3","all315_maf0.05","dedup188_mac2"))
  add(paste0(s2, " (length-matched power gate)"), s2,
      "allele count threshold, gate held fixed", sub("all315_|dedup188_", "", s2),
      power_mode = "cds_bp")
add("all factors loosened", "all315_mac1", "combination",
    "callab 200 + loose + exp_tot>=10 + pc1", callab_frac = 200/315, cutoff = "loose",
    power_min = 10, shrink = "pc1")

cat("=== running", length(grid), "parameter settings ===\n")
datasets <- unique(sapply(grid, function(g) g$data))
dsl <- setNames(lapply(datasets, mk), datasets)

res_list <- lapply(grid, function(g) {
  dd <- dsl[[g$data]]
  r <- run_funnel(dd$d, dd$n, callab_frac = g$callab_frac, shrink = g$shrink,
                  cutoff = g$cutoff, power_min = g$power_min, power_mode = g$power_mode,
                  dens_ref = DENS_REF, eval_rule = g$evaluable)
  m <- metrics(r, g$label, g$factor, g$level)
  m[, data_setting := g$data]
  m[, n_genomes := dd$n]
  m[, variant_density := r$dens]
  m[, median_tot_seg_evaluable := median(r$all[callable & evaluable, syn_seg + nonsyn_seg])]
  list(m = m, funnel = r$funnel, ranked = r$ranked)
})
M <- rbindlist(lapply(res_list, `[[`, "m"))
setcolorder(M, c("setting","factor","level","data_setting","n_genomes"))
print(M[, .(setting, rho = round(rho_genomewide, 3), tau = round(tau_genomewide, 3),
            rho_sl = round(rho_shortlist, 3), n_elig = n_eligible,
            J20 = round(jaccard_top20, 2), J50 = round(jaccard_top50, 2),
            J100 = round(jaccard_top100, 2), t1 = tier1_eligible)])
fwrite(M, file.path(OT, "rev_b2_settings_metrics.csv"))

FC <- rbindlist(lapply(seq_along(grid), function(i)
  data.table(setting = grid[[i]]$label, step = names(res_list[[i]]$funnel),
             n = as.integer(res_list[[i]]$funnel))))
fwrite(dcast(FC, setting ~ step, value.var = "n"), file.path(OT, "rev_b2_funnel_counts.csv"))

## ---------------------------------------------------------------- resampling
rs_names <- grep("^resample", unique(cnt$setting), value = TRUE)
cat("\n=== balanced lineage resampling:", length(rs_names), "draws ===\n")
cat("one genome drawn at random per lineage (188 lineages), published estimator/gates\n")

rs <- mclapply(rs_names, function(s) {
  dd <- mk(s)
  r <- run_funnel(dd$d, dd$n)
  m <- metrics(r, s, "lineage resample", "1 genome per lineage")
  pc <- r$all[locus_tag %in% shortlist, .(locus_tag, constraint_pctile, cstat,
                                          eligible = locus_tag %in% r$ranked$locus_tag)]
  pc[, draw := s]
  rk <- r$ranked[, .(locus_tag, rank)]
  rk[, draw := s]
  list(m = m, pc = pc, rk = rk)
}, mc.cores = NCORE)

RM <- rbindlist(lapply(rs, `[[`, "m"))
fwrite(RM, file.path(OT, "rev_b2_resample_metrics.csv"))
cat(sprintf("rho genome-wide: median %.3f  range %.3f-%.3f\n",
            median(RM$rho_genomewide), min(RM$rho_genomewide), max(RM$rho_genomewide)))
cat(sprintf("tau genome-wide: median %.3f  range %.3f-%.3f\n",
            median(RM$tau_genomewide), min(RM$tau_genomewide), max(RM$tau_genomewide)))
cat(sprintf("rho shortlist  : median %.3f  range %.3f-%.3f\n",
            median(RM$rho_shortlist), min(RM$rho_shortlist), max(RM$rho_shortlist)))
for (k in c(20, 50, 100)) {
  v <- RM[[paste0("jaccard_top", k)]]
  cat(sprintf("Jaccard top-%-3d: median %.3f  range %.3f-%.3f\n", k, median(v), min(v), max(v)))
}
cat(sprintf("n eligible     : median %.0f  range %.0f-%.0f\n",
            median(RM$n_eligible), min(RM$n_eligible), max(RM$n_eligible)))

PC <- rbindlist(lapply(rs, `[[`, "pc"))
RK <- rbindlist(lapply(rs, `[[`, "rk"))

cand_sum <- PC[, .(n_draws = .N,
                   pctile_median = median(constraint_pctile, na.rm = TRUE),
                   pctile_q05 = quantile(constraint_pctile, 0.05, na.rm = TRUE),
                   pctile_q95 = quantile(constraint_pctile, 0.95, na.rm = TRUE),
                   pctile_min = min(constraint_pctile, na.rm = TRUE),
                   pctile_max = max(constraint_pctile, na.rm = TRUE),
                   eligible_freq = mean(eligible)), by = locus_tag]
rk_sum <- RK[, .(rank_median = as.numeric(median(as.numeric(rank))),
                 rank_q05 = as.numeric(quantile(as.numeric(rank), .05)),
                 rank_q95 = as.numeric(quantile(as.numeric(rank), .95)),
                 n_draws_ranked = .N), by = locus_tag]
cand_sum <- merge(cand_sum, rk_sum, by = "locus_tag", all.x = TRUE)
cand_sum <- merge(pub_cand[, .(locus_tag, published_rank = rank, product,
                               published_pctile = constraint_pctile,
                               published_piNpiS = piNpiS_shrunk)],
                  cand_sum, by = "locus_tag", all.x = TRUE)
cand_sum[, tier1 := locus_tag %in% t1]
setorder(cand_sum, published_rank)
fwrite(cand_sum, file.path(OT, "rev_b2_candidate_percentiles.csv"))

t1tab <- cand_sum[tier1 == TRUE, .(locus_tag, product, published_rank, published_pctile,
                                   published_piNpiS, pctile_median, pctile_q05, pctile_q95,
                                   rank_median, rank_q05, rank_q95,
                                   eligible_freq, n_draws_ranked)]
setorder(t1tab, published_rank)
cat("\n=== Tier-1 retention across", length(rs_names), "balanced lineage resamples ===\n")
print(t1tab[, .(locus_tag, pub_rank = published_rank, pub_pct = published_pctile,
                pct_med = round(pctile_median, 1),
                pct_90CI = paste0(round(pctile_q05, 1), "-", round(pctile_q95, 1)),
                rank_med = rank_median,
                rank_90CI = paste0(rank_q05, "-", rank_q95),
                elig_freq = round(eligible_freq, 3))])
fwrite(t1tab, file.path(OT, "rev_b2_tier1_retention.csv"))

## Tier-1 eligibility across the deterministic parameter grid, per candidate
t1grid <- rbindlist(lapply(seq_along(grid), function(i) {
  rk <- res_list[[i]]$ranked
  data.table(setting = grid[[i]]$label, factor = grid[[i]]$factor, locus_tag = t1,
             eligible = t1 %in% rk$locus_tag,
             rank = rk$rank[match(t1, rk$locus_tag)])
}))
fwrite(dcast(t1grid, locus_tag ~ setting, value.var = "eligible"),
       file.path(OT, "rev_b2_tier1_eligibility_by_setting.csv"))
cat("\nTier-1 eligibility frequency across the", length(grid), "deterministic settings:\n")
print(t1grid[, .(elig_frac = round(mean(eligible), 3),
                 median_rank = as.numeric(median(as.numeric(rank), na.rm = TRUE))),
             by = locus_tag][order(-elig_frac)])

saveRDS(list(M = M, RM = RM, PC = PC, RK = RK, cand_sum = cand_sum, t1grid = t1grid,
             pub_cand = pub_cand, tier1 = tier1),
        file.path(Sys.getenv("REV_B2_SCRATCH",
          "/tmp/claude-1000/-mnt-HuangLab-26T-HuangLab-Data-Tar-Nine-Ta/459df96a-3261-4f14-8f94-382b2767e924/scratchpad"),
          "rev_b2_asp_results.rds"))
cat("\ndone\n")
