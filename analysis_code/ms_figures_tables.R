# Manuscript figures (PDF -> output/figure/ms_fig*.pdf; PNG copies for the Word draft)
# and manuscript tables (CSV -> output/table/ms_table*.csv).
# Inputs are existing project tables only; no analysis is re-run here.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel); library(scales)
})
# ggrepel places labels stochastically, so without a seed two runs of this script produce
# figures that differ in label placement. Fixed here so the deposited code redraws them exactly.
set.seed(20261012)
base   <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
td     <- file.path(base, "doc/Aspergillus/target_discovery_2026-08/tables")
figdir <- file.path(base, "output/figure")
tabdir <- file.path(base, "output/table")
pngdir <- file.path(base, "doc/manuscript_2026-10/figures")
dir.create(pngdir, showWarnings = FALSE, recursive = TRUE)

theme_ms <- theme_classic(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 11),
        axis.text = element_text(colour = "black"),
        legend.key.size = unit(3, "mm"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
col_good <- "#1b7a5a"; col_mid <- "#b07d12"; col_bad <- "#b2182b"; col_grey <- "grey70"

save_fig <- function(p, name, w, h) {
  ggsave(file.path(figdir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(pngdir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}

# ---------------------------------------------------------------- Figure 1: workflow
box <- function(x, y, w, h, label, fill) data.table(x, y, w, h, label, fill)
bx <- rbind(
  box(2.5, 9.0, 4.4, 1.0, "NCBI: 1,097 Aspergillus assemblies (11 species)\nQC: 934 pass (depth, N50, size, contigs)", "#e8eef4"),
  box(0.95, 7.2, 1.8, 1.2, "Phenotype-linked route\n55 A. fumigatus genomes\nwith numeric MIC\n(literature-mined)", "#f6e8e8"),
  box(0.95, 5.4, 1.8, 1.4, "LMM-GWAS (unitigs, SNPs)\n+ cyp51A positive control\n+ 31-gene panel\n=> underpowered", "#f6e8e8"),
  box(3.4, 7.2, 2.3, 1.2, "Phenotype-free route\n315 A. fumigatus genomes\n477,735 variant alleles\nper-gene piN/piS (9,823 genes)", "#e6f2ec"),
  box(3.4, 5.4, 2.3, 1.4, "Selectivity & tractability\nno human homolog, no paralog,\npan-Aspergillus, powered\n=> 815 candidates", "#e6f2ec"),
  box(3.4, 3.5, 2.3, 1.3, "Structure & chemistry (top 200)\nAlphaFold pLDDT, fpocket,\nPDB ligand-bound homologs\n=> 18 Tier-1", "#e6f2ec"),
  box(6.1, 7.2, 2.3, 1.2, "Second genus\n424 C. tropicalis isolates\nvs MYA-3404 (reads, diploid)\nper-gene piN/piS (6,254 genes)", "#e8e3f2"),
  box(6.1, 5.4, 2.3, 1.4, "Cross-genus intersection\n3,295 orthologue pairs\nselectivity vs spectrum\n=> 18 shared candidates", "#e8e3f2"),
  box(6.1, 3.5, 2.3, 1.3, "External benchmark\n25,536 C. albicans\nessentiality calls\n=> AUC 0.57-0.70", "#e8e3f2"),
  box(4.75, 2.05, 5.0, 0.85, paste0("Independent validation layer (added at revision)\n",
      "baseline-model benchmarking | ranking-stability resampling\n",
      "orthogroup phylogeny | pathway enrichment | power simulation"), "#f2eee4"),
  box(4.75, 0.78, 5.0, 0.85, "Evidence-graded candidate hypotheses, with their validation status stated", "#ece8f4")
)
ar <- data.table(x    = c(1.6, 3.2, 0.95, 3.4, 3.4, 3.4,  6.1, 6.1, 6.1,  4.75),
                 y    = c(8.5, 8.5, 6.6,  6.6, 4.7, 2.85, 6.6, 4.7, 2.85, 1.625),
                 xend = c(1.1, 3.4, 0.95, 3.4, 3.4, 3.4,  6.1, 6.1, 6.1,  4.75),
                 yend = c(7.8, 7.8, 6.1,  6.1, 4.15, 2.53, 6.1, 4.15, 2.53, 1.26))
f1 <- ggplot() +
  geom_rect(data = bx, aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = I(fill)),
            colour = "grey30", linewidth = 0.3) +
  geom_text(data = bx, aes(x, y, label = label), size = 2.4, lineheight = 0.95) +
  geom_segment(data = ar, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = 0.35) +
  annotate("segment", x = 1.85, xend = 2.2, y = 5.4, yend = 5.4, linetype = "dashed",
           arrow = arrow(length = unit(1.6, "mm")), linewidth = 0.3) +
  annotate("text", x = 2.02, y = 5.62, label = "motivates", size = 2.1, fontface = "italic") +
  annotate("segment", x = 4.6, xend = 4.9, y = 7.2, yend = 7.2, linetype = "dashed", linewidth = 0.3) +
  annotate("segment", x = 4.6, xend = 4.9, y = 5.4, yend = 5.4, linetype = "dashed", linewidth = 0.3) +
  coord_cartesian(xlim = c(-0.1, 7.4), ylim = c(0.2, 9.6), expand = FALSE) + theme_void()
save_fig(f1, "ms_fig1_workflow", 8.0, 6.0)

# ---------------------------------------------------------------- Figure 2: phenotype route
v <- fread(file.path(tabdir, "vrc_cyp51A_genotypes.csv"))
v[, cls := fifelse(TR == "TR46", "TR46/Y121F/T289A", fifelse(TR == "TR34", "TR34/L98H",
             fifelse(hotspot != "none", "coding hotspot only", "no cyp51A mechanism")))]
v[, cls := factor(cls, levels = c("no cyp51A mechanism", "coding hotspot only", "TR34/L98H", "TR46/Y121F/T289A"))]
short <- c("no cyp51A mechanism" = "None", "coding hotspot only" = "Hotspot\nonly", "TR34/L98H" = "TR34/\nL98H",
           "TR46/Y121F/T289A" = "TR46/\nY121F/T289A")
n_lab <- v[, .N, by = cls][, lab := paste0(short[as.character(cls)], "\n(n = ", N, ")")]
v <- merge(v, n_lab[, .(cls, lab)], by = "cls")
v[, lab := factor(lab, levels = n_lab[order(cls), lab])]
f2a <- ggplot(v, aes(lab, mic_used)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey40") +
  geom_jitter(aes(colour = R_call), width = 0.15, height = 0.06, size = 1.3, alpha = 0.85) +
  scale_y_continuous(trans = "log2", breaks = 2^(-2:6), labels = label_number(accuracy = 0.01, drop0trailing = TRUE)) +
  scale_colour_manual(values = c(R = col_bad, S = "grey45"), name = "EUCAST call") +
  labs(x = "cyp51A genotype class", y = "Voriconazole MIC (mg/L)") + theme_ms +
  theme(axis.text.x = element_text(size = 7.5), legend.position = c(0.2, 0.85))
s <- fread(file.path(tabdir, "vrc_gwas_scan_summary.csv"))
# Positive controls are READ from the association output, not hard-coded. The earlier version of
# this panel carried p = 3.3e-5 and 4.2e-4, which the revision withdrew: 3.3e-5 is the minimum P
# of the 46-genome binary scan and not a single-hypothesis test, and 4.2e-4 is in no output. Each
# control is now taken from the cohort that matches its phenotype's own genome-wide scan, so every
# point in this panel is comparable with a threshold drawn in the same panel: binary resistance
# from the 55-genome scan, log2 MIC from the 46-genome lineage-de-duplicated scan.
pc_file <- function(ph, coh) file.path(base, "analysis/vrc_gwas/04_association",
                                       sprintf("known_mech_%s_%s.tsv", ph, coh))
pc <- rbindlist(list(
  fread(pc_file("binary",  "all55"))[,   .(variant, p = `lrt-pvalue`, arm = "R/S, n = 55")],
  fread(pc_file("log2mic", "dedup46"))[, .(variant, p = `lrt-pvalue`, arm = "log2 MIC, n = 46")]))
pc <- pc[variant %in% c("cyp51A_TR46", "cyp51A_known_mech")]
pc[, label := paste0(fifelse(variant == "cyp51A_TR46", "cyp51A TR46", "any cyp51A mechanism"),
                     " (single test, ", arm, ")")]
pp <- rbind(
  s[, .(label, p = min_p, thr = threshold, type = "Genome-wide scan (best p)")],
  pc[, .(label, p, thr = NA_real_, type = "Known mechanism, same LMM")])
pp[, label := factor(label, levels = rev(label))]
f2b <- ggplot(pp, aes(-log10(p), label, colour = type)) +
  geom_point(data = pp[!is.na(thr)], aes(x = -log10(thr), y = label), shape = "|", size = 5, colour = "black", inherit.aes = FALSE) +
  geom_vline(xintercept = -log10(1.84e-7), linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("Genome-wide scan (best p)" = "grey35", "Known mechanism, same LMM" = col_bad), name = NULL) +
  scale_x_continuous(limits = c(0, 7.5)) +
  labs(x = expression(-log[10]~italic(p)), y = NULL) + theme_ms + theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2))
f2 <- (f2a | f2b) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A")
save_fig(f2, "ms_fig2_phenotype_route", 7.2, 3.4)

# ---------------------------------------------------------------- Figure 3: constraint metric
g <- fread(file.path(td, "constraint_all_genes.csv"))
known <- data.table(locus_tag = c("AFUA_5G04170","AFUA_6G12400","AFUA_5G09360","AFUA_4G04820","AFUA_6G05300",
                                  "AFUA_2G03700","AFUA_2G11010","AFUA_4G06890","AFUA_1G14870"),
                    name = c("hsp90","fks1","cnaA","erg25B","hapE","hmg1","DHODH (olorofim)","cyp51A (azoles)","gwt1 (fosmanogepix)"))
known <- merge(known, g[, .(locus_tag, piNpiS_shrunk)], by = "locus_tag")
known[, grp := fifelse(piNpiS_shrunk < 0.12, "a", fifelse(piNpiS_shrunk < 0.4, "b", "c"))]
gg <- g[!is.na(piNpiS_shrunk)]
gg[, x := pmin(piNpiS_shrunk, 1.5)]
med <- median(gg$piNpiS_shrunk)
f3a <- ggplot(gg, aes(x)) + geom_histogram(binwidth = 0.025, fill = "grey85", colour = NA, boundary = 0) +
  geom_vline(xintercept = med, linetype = "dashed", colour = "grey40") +
  geom_point(data = known, aes(x = piNpiS_shrunk, y = 0, colour = grp), size = 1.8, inherit.aes = FALSE) +
  geom_text_repel(data = known, aes(x = piNpiS_shrunk, y = 0, label = name, colour = grp), inherit.aes = FALSE,
                  nudge_y = 260, direction = "x", angle = 90, hjust = 0, size = 2.4, segment.size = 0.2,
                  segment.colour = "grey50", max.overlaps = Inf, min.segment.length = 0) +
  scale_colour_manual(values = c(a = col_good, b = col_mid, c = col_bad), guide = "none") +
  labs(x = expression(pi[N]/pi[S]~"(shrunk; values > 1.5 binned at 1.5)"), y = "Genes") + theme_ms
ev <- gg[mean_cov_genomes >= 250]
ev[, bin := cut(rank(piNpiS_shrunk, ties.method = "first"), 20, labels = FALSE)]
hb <- ev[, .(pct = 100 * mean(has_human_homolog)), by = bin]
f3b <- ggplot(hb, aes(bin, pct)) + geom_col(fill = col_mid, width = 0.8) +
  geom_hline(yintercept = 100 * mean(ev$has_human_homolog), linetype = "dashed") +
  labs(x = "Constraint bin (1 = most constrained, 20 = least)", y = "Genes with human homolog (%)") + theme_ms
k <- fread(file.path(td, "kinome_vanrhijn2024.csv"))
k[, locus_tag := toupper(gsub("\\s", "", locus_tag))]   # same filter as td_figures.R (Fig 5 of the report)
k <- merge(k, g[, .(locus_tag, piNpiS_shrunk, mean_cov_genomes, tot_seg)], by = "locus_tag")[mean_cov_genomes >= 250 & tot_seg >= 3]
k[, ess := factor(fifelse(homokaryotic_null_generated == "No", "No null obtained\n(essential)", "Null obtained\n(non-essential)"),
                  levels = c("No null obtained\n(essential)", "Null obtained\n(non-essential)"))]
wt <- wilcox.test(piNpiS_shrunk ~ ess, data = k)
auc <- with(k, { r <- rank(-piNpiS_shrunk); n1 <- sum(ess == levels(ess)[1]); n0 <- sum(ess != levels(ess)[1])
  (sum(r[ess == levels(ess)[1]]) - n1 * (n1 + 1) / 2) / (n1 * n0) })
f3c <- ggplot(k, aes(ess, piNpiS_shrunk)) + geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.6, colour = col_good) +
  annotate("text", x = 1.5, y = max(k$piNpiS_shrunk) * 1.02, size = 2.5,
           label = sprintf("n = %d vs %d; Wilcoxon p = %.1e; AUC = %.3f",
                           sum(k$ess == levels(k$ess)[1]), sum(k$ess != levels(k$ess)[1]), wt$p.value, auc)) +
  labs(x = NULL, y = expression(pi[N]/pi[S])) + theme_ms
dd <- fread(file.path(td, "constraint_dedup_comparison.csv"))[!is.na(all315) & !is.na(dedup)]
rho <- cor(dd$all315, dd$dedup, method = "spearman")
f3d <- ggplot(dd, aes(pmin(all315, 2), pmin(dedup, 2))) + geom_point(size = 0.3, alpha = 0.25) +
  geom_abline(colour = col_bad, linewidth = 0.4) +
  annotate("text", x = 0.05, y = 1.9, hjust = 0, size = 2.6, label = sprintf("Spearman rho = %.3f", rho)) +
  labs(x = expression(pi[N]/pi[S]~"all 315 genomes"), y = expression(pi[N]/pi[S]~"188 lineage representatives")) +
  coord_equal() + theme_ms
f3 <- (f3a | f3b) / (f3c | f3d) + plot_annotation(tag_levels = "A")
save_fig(f3, "ms_fig3_constraint_validation", 7.2, 6.4)
cat(sprintf("Fig3: median %.3f; human-homolog bin1 %.1f%% bin20 %.1f%% baseline %.1f%%; kinase AUC %.3f p %.2e; rho %.3f\n",
            med, hb[bin == 1, pct], hb[bin == 20, pct], 100 * mean(ev$has_human_homolog), auc, wt$p.value, rho))

# ---------------------------------------------------------------- Figure 4: prioritisation
# Funnel counts are READ, never typed. rev_c1_canonical_filter_audit.R derives every count
# gene by gene from the canonical inputs and writes this table; advisor round 2 (WP1) requires
# one reproducible table behind every published count.
fc  <- fread(file.path(tabdir, "rev_c1_funnel_counts_canonical.csv"))
SHORT <- c("(none)"          = "Protein-coding genes",
           "c1_callable"     = "Callable (>=250/315 genomes)",
           "c2_evaluable"    = "Evaluable",
           "c3_constrained"  = "Constrained",
           "c4_no_human_sen" = "No human homolog",
           "c5_no_paralog"   = "No A. fumigatus paralog",
           "c6_conserved"    = "Conserved in >=8/9 Aspergillus",
           "c7_not_known"    = "Not an already-drugged target",
           "c8_len_ok"       = "Protein >=100 aa",
           "c9_powered"      = "Powered (>=15 expected variants)")
fun <- data.table(step = unname(SHORT[fc$criterion]), n = fc$n_published_rule)
# the same funnel re-run end to end with the high-sensitivity human screen substituted
fun <- rbind(fun, data.table(step = "Human re-screen at high sensitivity",
                             n = fc[.N, n_corrected_rule]))
stopifnot(!any(is.na(fun$step)), nrow(fun) == 11L)
fun[, step := factor(step, levels = rev(step))]
f4a <- ggplot(fun, aes(n, step)) + geom_col(fill = c(rep("grey60", nrow(fun) - 1L), col_bad), width = 0.7) +
  geom_text(aes(label = comma(n)), hjust = -0.1, size = 2.5) +
  scale_x_continuous(limits = c(0, 11500), labels = comma) + labs(x = "Genes remaining", y = NULL) + theme_ms
ef <- data.table(step = c("Top-200 candidates", "AlphaFold model", "Well folded (pLDDT >= 70)",
                          "Druggable pocket (>= 0.5)", "Robust to de-duplication", "Ligand-bound PDB homolog (Tier 1)"),
                 n = c(200, 190, 113, 77, 69, 18))
ef[, step := factor(step, levels = rev(step))]
f4b <- ggplot(ef, aes(n, step)) + geom_col(fill = col_good, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.5) + scale_x_continuous(limits = c(0, 230)) +
  labs(x = "Candidates remaining", y = NULL) + theme_ms
fa <- fread(file.path(td, "candidates_full_assessment.csv"))[!is.na(mean_plddt)]
t1 <- fread(file.path(td, "tier1_shortlist.csv"))
fa[, tier := fifelse(locus_tag %in% t1$locus_tag, "Tier 1", "Other top-200")]
lab_ids <- c(AFUA_2G01220 = "GTP cyclohydrolase II", AFUA_3G10310 = "acetolactate synthase",
             AFUA_6G03540 = "malate synthase", AFUA_5G09750 = "NCS2 permease", AFUA_8G05210 = "SAM-MTase (rank 1)")
fa[, lab := lab_ids[locus_tag]]
f4c <- ggplot(fa, aes(mean_plddt, pmax(best_druggability, 0, na.rm = TRUE))) +
  annotate("rect", xmin = 70, xmax = 100, ymin = 0.5, ymax = 1.02, fill = "#e6f2ec") +
  geom_point(aes(colour = tier, size = tier), alpha = 0.8) +
  geom_text_repel(aes(label = lab), size = 2.4, min.segment.length = 0, box.padding = 0.4, na.rm = TRUE, max.overlaps = Inf) +
  scale_colour_manual(values = c("Tier 1" = col_good, "Other top-200" = "grey65"), name = NULL) +
  scale_size_manual(values = c("Tier 1" = 1.8, "Other top-200" = 0.9), name = NULL) +
  labs(x = "Mean AlphaFold pLDDT", y = "Best fpocket druggability score") + theme_ms +
  theme(legend.position = c(0.2, 0.9))
f4 <- (f4a / f4b + plot_layout(heights = c(1.6, 1))) | f4c
f4 <- f4 + plot_layout(widths = c(1, 1.05)) + plot_annotation(tag_levels = "A")
save_fig(f4, "ms_fig4_prioritisation", 7.4, 4.6)

# ---------------------------------------------------------------- Figure 5: proteomics
# DEFERRED to a follow-up paper (advisor's decision, 2026-10-11): the in-house LC-MS/MS
# comparison has one culture per condition and no biological replication, so it is held
# back until a properly replicated dataset is available. The code is kept intact; set
# MAKE_PROTEOMICS_FIG <- TRUE to regenerate the volcano plot for the follow-up analysis.
MAKE_PROTEOMICS_FIG <- FALSE
if (MAKE_PROTEOMICS_FIG) {
pr <- fread(file.path(base, "analysis/proteomics/proteomics_with_constraint.csv"))[!is.na(log2ratio) & !is.na(adj_p)]
pr <- pr[, .SD[which.min(adj_p)], by = locus_tag]
pr[, capped := abs(log2ratio) >= 6.6]
pr[, sig := adj_p < 0.05 & abs(log2ratio) >= 1]
pr[, cls := fifelse(capped, "on/off (capped ratio)", fifelse(sig & log2ratio > 0, "induced", fifelse(sig, "repressed", "n.s.")))]
hl <- c(AFUA_2G05060 = "alternative oxidase (AOX)", AFUA_4G06890 = "cyp51A", AFUA_7G04500 = "ATP-PRT (HisG)",
        AFUA_1G13740 = "AROM", AFUA_5G13330 = "mL67", AFUA_5G04170 = "hsp90")
pr[, lab := hl[locus_tag]]
f5 <- ggplot(pr, aes(log2ratio, -log10(adj_p))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
  geom_point(aes(colour = cls, shape = capped), size = 1.2, alpha = 0.8) +
  geom_text_repel(aes(label = lab), size = 2.5, na.rm = TRUE, min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf) +
  scale_colour_manual(values = c(induced = col_bad, repressed = "#2166ac", "n.s." = "grey70", "on/off (capped ratio)" = col_mid), name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  labs(x = expression(log[2]~"(voriconazole / untreated)"), y = expression(-log[10]~"adjusted"~italic(p))) + theme_ms +
  theme(legend.position = c(0.15, 0.85))
save_fig(f5, "ms_fig5_proteomics", 5.2, 4.0)
cat(sprintf("Fig5: %d genes; induced %d; repressed %d; capped %d\n", nrow(pr), pr[cls=="induced",.N], pr[cls=="repressed",.N], pr[capped==TRUE,.N]))
}

# ---------------------------------------------------------------- Tables (CSV)
t2 <- merge(known[, .(locus_tag, name)], g[, .(locus_tag, piNpiS_shrunk, has_human_homolog, has_paralog, n_species)], by = "locus_tag")
t2[, pctile := round(100 * ecdf(gg$piNpiS_shrunk)(piNpiS_shrunk) , 1)]
t2[, pctile := round(100 - pctile, 1)]  # percent of genes LESS constrained
setorder(t2, piNpiS_shrunk)
fwrite(t2, file.path(tabdir, "ms_table2_known_targets.csv"))
t3 <- t1[, .(rank, locus_tag, product, piNpiS = round(piNpiS_shrunk, 3), constraint_pctile,
             syn_seg, nonsyn_seg, plddt = round(mean_plddt, 1), pocket = best_druggability,
             pdb = best_pdb_id, pdb_identity = best_pdb, ligand = top_ligand, human_homolog_sensitive = human_sens)]
setorder(t3, rank)
fwrite(t3, file.path(tabdir, "ms_table3_tier1.csv"))
print(t2); print(t3)

# ---------------------------------------------------------------- combined main-text figure
# Main Figure 1 for the concise version: study design beside the filtering funnel.
fM1 <- (f1 / f4a) + plot_layout(heights = c(1.55, 1)) + plot_annotation(tag_levels = "A")
save_fig(fM1, "ms_figM1_design_funnel", 7.3, 7.6)

# ---------------------------------------------------------------- ASM/mSystems figure files
# ASM asks for each figure as one file, RGB, 300-600 dpi, <= 7 x 9 in, Arial/Helvetica/Times.
asmdir <- file.path(base, "doc/manuscript_2026-10/figures_mSystems")
dir.create(asmdir, showWarnings = FALSE, recursive = TRUE)
asm_tiff <- function(p, n, w, h) ggsave(file.path(asmdir, paste0(n, ".tiff")), p,
                                        width = min(w, 7), height = min(h, 9), dpi = 300,
                                        compression = "lzw", bg = "white")
asm_tiff(fM1, "FIG1", 7.0, 7.3)
asm_tiff(f3,  "FIG2", 7.0, 6.4)
if (MAKE_PROTEOMICS_FIG) asm_tiff(f5, "FIG5", 5.2, 4.0)
asm_tiff(f2,  "FIG_S1", 7.0, 3.4)
asm_tiff(f4,  "FIG_S2", 7.0, 4.6)
