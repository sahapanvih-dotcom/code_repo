# Main-text Figure 5 for the revised manuscript: what the added checks returned.
# (A) per-dataset ROC-AUC with 95% CI for the four prespecified models
# (B) paired delta-AUC of constraint against each baseline, and of the integrated
#     model against covariates, on identical gene-family-blocked resamples
# (C) ranking stability: top-20 Jaccard overlap against the published ranking,
#     by parameter factor, plus the lineage-resampling distribution
# Inputs: output/table/rev_b1_performance.csv, rev_b1_paired_deltas.csv,
#         output/table/rev_b2_settings_metrics.csv, rev_b2_resample_metrics.csv
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
base   <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
tabdir <- file.path(base, "output/table")
figdir <- file.path(base, "output/figure")
pngdir <- file.path(base, "doc/manuscript_2026-10/figures")
asmdir <- file.path(base, "doc/manuscript_2026-10/figures_mSystems")

theme_ms <- theme_classic(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 11),
        axis.text = element_text(colour = "black"),
        legend.key.size = unit(3, "mm"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))

perf <- fread(file.path(tabdir, "rev_b1_performance.csv"))[rule == "strict"]
perf <- perf[model %in% c("M1_constraint","M2_similarity","M3_covariates","M4_integrated")]
mlab <- c(M1_constraint = "constraint", M2_similarity = "similarity",
          M3_covariates = "covariates", M4_integrated = "integrated")
perf[, model_lab := factor(mlab[model], levels = unname(mlab))]
perf[, src := gsub("_", " ", source)]
# drop the two n-limited GRACEv3 panels from the figure but keep them in the table
perf_f <- perf[!grepl("GRACEv3", source)]

pA <- ggplot(perf_f, aes(src, roc_auc, colour = model_lab)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_pointrange(aes(ymin = roc_lo, ymax = roc_hi),
                  position = position_dodge(width = 0.62), size = 0.22, linewidth = 0.4) +
  facet_wrap(~ species_lab, nrow = 1) +
  scale_colour_manual(values = c(constraint = "#1b7a5a", similarity = "#b07d12",
                                 covariates = "#4575b4", integrated = "#b2182b"), name = NULL) +
  labs(x = NULL, y = "ROC-AUC (95% CI)") + theme_ms +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

dl <- fread(file.path(tabdir, "rev_b1_paired_deltas.csv"))[rule == "strict"]
keep <- c("M1_constraint - B_random", "M1_constraint - M2_similarity",
          "M1_constraint - M3_covariates", "M4_integrated - M3_covariates")
klab <- c("M1_constraint - B_random"        = "vs random",
          "M1_constraint - M2_similarity"   = "vs similarity",
          "M1_constraint - M3_covariates"   = "vs covariates",
          "M4_integrated - M3_covariates"   = "integrated model\nvs covariates")
dl <- dl[contrast %in% keep][!grepl("GRACEv3", source)]
dl[, lab := factor(klab[contrast], levels = unname(klab))]
dl[, sig := d_roc_lo > 0 | d_roc_hi < 0]

pB <- ggplot(dl, aes(lab, d_roc)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_pointrange(aes(ymin = d_roc_lo, ymax = d_roc_hi, colour = sig),
                  position = position_jitter(width = 0.16, height = 0, seed = 1),
                  size = 0.2, linewidth = 0.35) +
  facet_wrap(~ species_lab, nrow = 1) +
  scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey65"),
                      labels = c(`TRUE` = "CI excludes 0", `FALSE` = "CI spans 0"), name = NULL) +
  labs(x = "constraint compared with:", y = expression(paired~Delta*"AUC")) + theme_ms +
  theme(legend.position = "top", axis.text.x = element_text(angle = 20, hjust = 1))

st <- fread(file.path(tabdir, "rev_b2_settings_metrics.csv"))
st <- st[setting != "baseline (published)"]
fl <- c("shrinkage estimator" = "shrinkage", "expected-variant gate" = "power gate",
        "minor-allele filter" = "allele filter", "callability" = "callability",
        "constrained-gene cutoff" = "cutoff", "lineage sampling" = "lineage sampling",
        "combination" = "combination", "evaluable gate" = "evaluable gate")
st[, fac := ifelse(factor %in% names(fl), fl[factor], factor)]
st[, fac := gsub(", gate held fixed", " (gate fixed)", fac)]
rs <- fread(file.path(tabdir, "rev_b2_resample_metrics.csv"))
rs_dt <- data.table(fac = "lineage resampling", jaccard_top20 = rs$jaccard_top20)
st_dt <- st[, .(fac, jaccard_top20)]
both <- rbind(st_dt, rs_dt)
ord  <- both[, .(m = median(jaccard_top20, na.rm = TRUE)), by = fac][order(m)]
both[, fac := factor(fac, levels = ord$fac)]

pC <- ggplot(both, aes(fac, jaccard_top20)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, fill = "grey92", width = 0.6) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = NULL, y = "top-20 Jaccard overlap with the published ranking") + theme_ms

f5 <- (pA / pB / pC) + plot_layout(heights = c(1.15, 1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(figdir, "ms_fig5_benchmark_stability.pdf"), f5,
       width = 7.0, height = 8.2, device = cairo_pdf)
ggsave(file.path(pngdir, "ms_fig5_benchmark_stability.png"), f5,
       width = 7.0, height = 8.2, dpi = 300, bg = "white")
ggsave(file.path(asmdir, "FIG5.tiff"), f5, width = 7.0, height = 8.2,
       dpi = 300, compression = "lzw", bg = "white")
cat("Figure 5 written: PDF, PNG and FIG5.tiff\n")
cat("panel A datasets:", length(unique(perf_f$source)), "| panel B points:", nrow(dl),
    "| panel C groups:", length(unique(both$fac)), "\n")
