#!/usr/bin/env Rscript
# rev_b2: rank-stability figures + two diagnostics a reviewer would ask for.
# Outputs: output/figure/rev_b2_rank_stability.pdf
#          output/figure/rev_b2_crossgenus_matched.pdf
#          output/table/rev_b2_top_of_list_diagnostics.csv
#          output/table/rev_b2_setlevel_jaccard.csv
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
TD <- "analysis/target_discovery"; OT <- "output/table"; OF <- "output/figure"
SD <- Sys.getenv("REV_B2_SCRATCH",
  "/tmp/claude-1000/-mnt-HuangLab-26T-HuangLab-Data-Tar-Nine-Ta/459df96a-3261-4f14-8f94-382b2767e924/scratchpad")

A <- readRDS(file.path(SD, "rev_b2_asp_results.rds"))
X <- readRDS(file.path(SD, "rev_b2_crossgenus_results.rds"))
M <- A$M; RM <- A$RM; PC <- A$PC; cand <- A$cand_sum; t1 <- A$tier1$locus_tag

th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(colour = "black"),
        plot.tag = element_text(face = "bold", size = 10),
        legend.key.size = unit(3, "mm"), strip.background = element_blank(),
        legend.position = "bottom", legend.title = element_text(size = 7))
col_bad <- "#b2182b"; col_good <- "#1b7a5a"; col_mid <- "#2166ac"

## -------------------------------------------------- panel A/B: settings grid
S <- copy(M)[setting != "baseline (published)"]
S[, fac := factor(factor, levels = c("allele count threshold",
                                     "allele count threshold, gate held fixed",
                                     "population sampling", "callability threshold",
                                     "shrinkage estimator", "constrained cutoff",
                                     "expected-variant gate", "evaluable gate", "combination"))]
setorder(S, fac, -rho_genomewide)
S[, setting := factor(setting, levels = unique(setting))]

pA <- ggplot(S, aes(rho_genomewide, setting, colour = fac)) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50") +
  geom_segment(aes(x = 1, xend = rho_genomewide, yend = setting), linewidth = 0.3) +
  geom_point(size = 1.4) +
  geom_point(aes(x = rho_shortlist), shape = 1, size = 1.4) +
  scale_x_continuous(limits = c(0.65, 1.02), breaks = seq(0.7, 1, 0.1)) +
  scale_colour_brewer(palette = "Set1", name = NULL, guide = guide_legend(ncol = 2)) +
  labs(x = expression("Spearman"~rho~"vs published ranking (filled = genome-wide, open = shortlist)"),
       y = NULL) + th + theme(axis.text.y = element_text(size = 5.6))

Jl <- melt(S[, .(setting, fac, `20` = jaccard_top20, `50` = jaccard_top50,
                 `100` = jaccard_top100)],
           id.vars = c("setting", "fac"), variable.name = "k", value.name = "J")
Jl[, k := factor(k, levels = c("20", "50", "100"))]
pB <- ggplot(Jl, aes(J, setting, colour = k, shape = k)) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50") +
  geom_point(size = 1.3, alpha = 0.85) +
  scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = c(`20` = col_bad, `50` = "#ef8a62", `100` = col_mid),
                      name = "top-k") +
  scale_shape_manual(values = c(16, 17, 15), name = "top-k") +
  labs(x = "Jaccard overlap with the published top-k", y = NULL) +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## -------------------------------------------------- panel C: lineage resampling
RL <- rbind(
  data.table(metric = "rho genome-wide", v = RM$rho_genomewide),
  data.table(metric = "rho shortlist", v = RM$rho_shortlist),
  data.table(metric = "Kendall tau genome-wide", v = RM$tau_genomewide),
  data.table(metric = "Jaccard top-20", v = RM$jaccard_top20),
  data.table(metric = "Jaccard top-50", v = RM$jaccard_top50),
  data.table(metric = "Jaccard top-100", v = RM$jaccard_top100))
RL[, metric := factor(metric, levels = unique(metric))]
pC <- ggplot(RL, aes(v, metric)) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50") +
  geom_boxplot(width = 0.55, outlier.size = 0.4, linewidth = 0.3,
               fill = "#9ecae1", colour = "grey25") +
  scale_x_continuous(limits = c(0.5, 1.02), breaks = seq(0.5, 1, 0.1)) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("%d balanced lineage resamples (1 genome per lineage, n = %d)",
                          nrow(RM), RM$n_genomes[1])) +
  th + theme(plot.subtitle = element_text(size = 6.8))

## -------------------------------------------------- panel D: Tier-1 percentile spread
T1 <- PC[locus_tag %in% t1]
ord <- cand[locus_tag %in% t1][order(published_rank)]
T1[, locus_tag := factor(locus_tag, levels = rev(ord$locus_tag))]
ordf <- copy(ord)[, locus_tag := factor(locus_tag, levels = rev(ord$locus_tag))]
pD <- ggplot(T1, aes(constraint_pctile, locus_tag)) +
  geom_boxplot(width = 0.55, outlier.size = 0.3, linewidth = 0.25,
               fill = "#c7e9c0", colour = "grey25") +
  geom_point(data = ordf, aes(published_pctile, locus_tag), colour = col_bad,
             shape = 124, size = 2.2) +
  labs(x = "constraint percentile across resamples (red bar = published value)", y = NULL,
       subtitle = "Tier-1 candidates; eligibility retained in 100/100 draws") +
  th + theme(axis.text.y = element_text(size = 5.6), plot.subtitle = element_text(size = 6.8))

## -------------------------------------------------- panel E: candidate rank spread
RK <- A$RK
rkj <- merge(RK, cand[, .(locus_tag, published_rank)], by = "locus_tag")
rkj <- rkj[published_rank <= 200]
pE <- ggplot(rkj, aes(published_rank, rank)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  geom_point(size = 0.25, alpha = 0.12, colour = "grey30") +
  geom_point(data = rkj[locus_tag %in% t1], size = 0.3, alpha = 0.3, colour = col_bad) +
  labs(x = "published rank", y = "rank in each lineage resample",
       subtitle = "published top-200 only; red = Tier-1") +
  th + theme(plot.subtitle = element_text(size = 6.8))

fig1 <- (pA | pB) / (pC | pD) / (pE | plot_spacer()) +
  plot_layout(heights = c(1.5, 1.1, 0.85)) +
  plot_annotation(tag_levels = "A",
    title = "rev_b2: stability of the A. fumigatus constraint ranking",
    theme = theme(plot.title = element_text(size = 10, face = "bold")))
ggsave(file.path(OF, "rev_b2_rank_stability.pdf"), fig1, width = 9.5, height = 11,
       device = cairo_pdf)
cat("wrote", file.path(OF, "rev_b2_rank_stability.pdf"), "\n")

## -------------------------------------------------- cross-genus figure
R <- copy(X$R)
R[, lab := factor(comparison, levels = rev(comparison))]
R[, grp := fifelse(grepl("^PUBLISHED", comparison), "published",
            fifelse(grepl("^MATCHED", comparison), "fully matched", "one axis matched"))]
pF <- ggplot(R, aes(spearman_rho, lab, colour = grp)) +
  geom_vline(xintercept = X$rho_pub, linetype = "dashed", colour = col_bad, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = rho_lo95, xmax = rho_hi95), height = 0, linewidth = 0.4) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(published = col_bad, `fully matched` = col_good,
                                 `one axis matched` = col_mid), name = NULL) +
  scale_x_continuous(limits = c(0.1, 0.36)) +
  labs(x = expression("cross-genus Spearman"~rho~"(95% bootstrap CI over orthologue pairs)"),
       y = NULL, subtitle = sprintf("dashed line = published rho = %.3f (n = %d pairs)",
                                    X$rho_pub, X$n_pub)) +
  th + theme(axis.text.y = element_text(size = 6), plot.subtitle = element_text(size = 6.8))

D <- X$D
pG <- ggplot(D, aes(factor(ct_n_iso), spearman_rho, colour = estimator)) +
  geom_hline(yintercept = X$rho_pub, linetype = "dashed", colour = col_bad, linewidth = 0.4) +
  geom_boxplot(width = 0.6, outlier.size = 0.4, linewidth = 0.3,
               position = position_dodge(0.7)) +
  scale_colour_manual(values = c(raw = col_mid, pc0.5 = col_good), name = "matched estimator") +
  labs(x = "C. tropicalis isolates retained (lineage representatives)",
       y = expression("cross-genus Spearman"~rho),
       subtitle = "sampling-depth matching: both genera deduplicated, estimator matched") +
  th + theme(plot.subtitle = element_text(size = 6.8))

REL <- melt(X$REL[, .(estimator, `A. fumigatus` = af_reliability_SB,
                      `C. tropicalis` = ct_reliability_SB)],
            id.vars = "estimator", variable.name = "genus", value.name = "rel")
pH <- ggplot(REL, aes(rel, genus, fill = estimator)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50") +
  scale_x_continuous(limits = c(0, 1.05)) +
  scale_fill_manual(values = c(raw = col_mid, pc0.5 = col_good), name = "estimator") +
  labs(x = "within-genus split-half reliability (Spearman-Brown)", y = NULL,
       subtitle = "estimator noise is NOT the reason rho is low") +
  th + theme(plot.subtitle = element_text(size = 6.8))

fig2 <- pF / (pG | pH) + plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = "A",
    title = "rev_b2: cross-genus constraint transfer under matched estimator and sampling",
    theme = theme(plot.title = element_text(size = 10, face = "bold")))
ggsave(file.path(OF, "rev_b2_crossgenus_matched.pdf"), fig2, width = 8.5, height = 7.5,
       device = cairo_pdf)
cat("wrote", file.path(OF, "rev_b2_crossgenus_matched.pdf"), "\n")

## -------------------------------------------------- diagnostic 1: top of the list
# How much of the head of the published ranking is driven by genes with (nearly) no
# observed nonsynonymous variants, where the shrunk ratio is set by the pseudocount?
pub_cand <- A$pub_cand
dg <- rbindlist(lapply(c(20, 50, 100, 200, 815), function(k) {
  h <- head(pub_cand, k)
  data.table(top_k = k,
             n = nrow(h),
             frac_nonsyn_seg_0 = mean(h$nonsyn_seg == 0),
             frac_nonsyn_seg_le1 = mean(h$nonsyn_seg <= 1),
             frac_nonsyn_seg_le2 = mean(h$nonsyn_seg <= 2),
             median_nonsyn_seg = median(h$nonsyn_seg),
             median_tot_seg = median(h$tot_seg),
             median_exp_tot = median(h$exp_tot))
}))
cat("\n=== diagnostic: observed nonsynonymous variants in the head of the published list ===\n")
print(dg)
fwrite(dg, file.path(OT, "rev_b2_top_of_list_diagnostics.csv"))

## -------------------------------------------------- diagnostic 2: set vs order
# Top-k Jaccard mixes two things: who is ELIGIBLE and how the eligible genes are
# ORDERED. The set-level Jaccard of the whole eligible list separates them.
setj <- M[, .(setting, factor, level, n_eligible,
              jaccard_eligible_set = round(jaccard_eligible_set, 3),
              frac_published_shortlist_retained = round(frac_published_shortlist_retained, 3),
              jaccard_top20 = round(jaccard_top20, 3),
              top20_median_nonsyn_seg)]
cat("\n=== eligible-SET overlap vs top-20 ORDER overlap ===\n")
print(setj)
fwrite(setj, file.path(OT, "rev_b2_setlevel_jaccard.csv"))

pI <- ggplot(M[setting != "baseline (published)"],
             aes(jaccard_eligible_set, jaccard_top20, colour = factor)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  geom_point(size = 1.5) +
  scale_colour_brewer(palette = "Set1", name = NULL,
                      guide = guide_legend(ncol = 2)) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Jaccard of the whole eligible shortlist (membership)",
       y = "Jaccard of the top 20 (order)",
       subtitle = "below the diagonal = the ORDER moves more than the MEMBERSHIP") + th +
  theme(plot.subtitle = element_text(size = 6.8))
ggsave(file.path(OF, "rev_b2_membership_vs_order.pdf"), pI, width = 5.2, height = 5,
       device = cairo_pdf)
cat("wrote", file.path(OF, "rev_b2_membership_vs_order.pdf"), "\n")
cat("\ndone\n")
