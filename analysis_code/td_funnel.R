#!/usr/bin/env Rscript
# Drug-target triage funnel for A. fumigatus.
# Constraint alone is not a target: the most constrained genes in any genome are the ones
# conserved across all life (ribosome, proteasome), which are exactly the ones the host shares.
# This applies the selectivity filters that turn constraint into a target list.
suppressPackageStartupMessages({library(data.table)})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
TD <- "analysis/target_discovery"

d <- fread(file.path(TD,"constraint_per_gene.tsv"))

## --- two complementary constraint statistics -------------------------------
# A: nonsynonymous depletion. Of the variants that DID occur in this gene, are too few
#    of them protein-changing? Internally controls for local mutation rate and gene length.
d[, tot_seg := syn_seg + nonsyn_seg]
d[, exp_nonsyn_frac := nonsyn_sites/(nonsyn_sites + syn_sites)]
d[, p_nonsyn_depl := ifelse(tot_seg >= 3,
      pbinom(nonsyn_seg, tot_seg, exp_nonsyn_frac), NA_real_)]
# shrunk ratio, defined even when syn_seg is small
d[, piNpiS_shrunk := ((nonsyn_seg + 0.5)/nonsyn_sites) / ((syn_seg + 0.5)/syn_sites)]

# B: total depletion. Given genome-wide variant density, does this gene carry far fewer
#    variants than its length predicts? Catches genes that barely vary at all.
cds_bp   <- sum(d$nonsyn_sites + d$syn_sites)
dens     <- sum(d$tot_seg)/cds_bp
d[, exp_tot := (nonsyn_sites + syn_sites)*dens]
d[, p_tot_depl := ppois(tot_seg, exp_tot)]

## --- evidence layers --------------------------------------------------------
hum <- fread("analysis/target_discovery/hom/af293_vs_human.tsv",
             col.names=c("protein_id","human_hit","pident","len","evalue","bits","qcov","scov"))
hum <- hum[, .SD[which.min(evalue)], by=protein_id]
par <- fread("analysis/target_discovery/hom/has_paralog.txt", header=FALSE)$V1
xs  <- fread("analysis/target_discovery/hom/crossspecies_besthits.tsv",
             header=FALSE, col.names=c("protein_id","species"))
xsn <- xs[, .(n_species=uniqueN(species)), by=protein_id]

d <- merge(d, hum[, .(protein_id, human_pident=pident, human_qcov=qcov, human_evalue=evalue)],
           by="protein_id", all.x=TRUE)
d <- merge(d, xsn, by="protein_id", all.x=TRUE)
d[is.na(n_species), n_species := 0L]
d[, has_paralog := protein_id %in% par]
d[, has_human_homolog := !is.na(human_evalue)]

## --- already-drugged / known targets, excluded from a NOVEL list ------------
known <- c("AFUA_4G06890","AFUA_7G03740","AFUA_6G12400","AFUA_2G11010","AFUA_1G14870",
           "AFUA_2G03700","AFUA_1G14330","AFUA_2G02690","AFUA_5G04170","AFUA_5G09360",
           "AFUA_4G04820","AFUA_6G05300","AFUA_2G01260","AFUA_3G09850")
d[, known_target := locus_tag %in% known]

## --- the funnel -------------------------------------------------------------
step <- function(dt, label) { cat(sprintf("%-52s %5d\n", label, nrow(dt))); dt }
cat("\n=== FUNNEL ===\n")
f0 <- step(d, "all Af293 protein-coding genes")
f1 <- step(f0[mean_cov_genomes >= 250], "callable: aligned in >=250 of 315 genomes")
f2 <- step(f1[tot_seg >= 3 | exp_tot >= 10],
           "evaluable: enough variants, or enough expected")
f3 <- step(f2[p_nonsyn_depl < 0.05 | p_tot_depl < 0.01],
           "constrained: nonsyn-depleted or variant-poor")
f4 <- step(f3[has_human_homolog == FALSE], "no human homolog (E>1e-3)  <- selectivity")
f5 <- step(f4[has_paralog == FALSE],       "no A. fumigatus paralog    <- no redundancy")
f6 <- step(f5[n_species >= 8],             "conserved in >=8 of 9 other Aspergillus")
f7 <- step(f6[known_target == FALSE],      "not an already-drugged target")
f8 <- step(f7[protein_len >= 100],         "protein >=100 aa (druggable size)")

## --- rank -------------------------------------------------------------------
# NOTE: -log10(p) scales with gene length, because long genes accumulate more variants and
# therefore more power. Ranking on the p-value alone floats long genes to the top regardless
# of how constrained they actually are. So: gate on POWER, then rank on EFFECT SIZE.
f8[, powered := exp_tot >= 15]                     # enough expected variants to say anything
f8[, constraint_score := -log10(pmin(
      ifelse(is.na(p_nonsyn_depl), 1, p_nonsyn_depl),
      ifelse(is.na(p_tot_depl),    1, p_tot_depl)) + 1e-300)]
# percentile of constraint among all callable genes, for interpretability
ref_dist <- d[mean_cov_genomes >= 250 & tot_seg >= 3]$piNpiS_shrunk
f8[, constraint_pctile := round(100*sapply(piNpiS_shrunk, function(x) mean(ref_dist > x, na.rm=TRUE)), 1)]

cat(sprintf("\n%-52s %5d\n", "adequately powered (>=15 expected variants)", sum(f8$powered)))
f9 <- f8[powered == TRUE]
setorder(f9, piNpiS_shrunk, -constraint_score)
f9[, rank := .I]

out <- f9[, .(rank, locus_tag, protein_id, protein_len, product,
              piN_piS, piNpiS_shrunk, constraint_pctile, syn_seg, nonsyn_seg, stop_seg,
              tot_seg, exp_tot=round(exp_tot,1), p_nonsyn_depl, p_tot_depl,
              constraint_score, n_species, mean_cov_genomes)]
fwrite(out, "doc/Aspergillus/target_discovery_2026-08/tables/candidate_targets_ranked.csv")
fwrite(d,   "doc/Aspergillus/target_discovery_2026-08/tables/constraint_all_genes.csv")

cat("\n=== TOP 25 CANDIDATES ===\n")
print(head(out[, .(rank, locus_tag, len=protein_len, piNpiS=round(piNpiS_shrunk,3),
                   pct=constraint_pctile, obs=tot_seg, exp=round(exp_tot,1), nsyn=nonsyn_seg,
                   product=substr(product,1,50))], 30))
cat("\nwrote candidate_targets_ranked.csv (", nrow(out), "candidates )\n")
