# Cross-genus intersection: A. fumigatus (315 assemblies, Af293 reference) vs
# C. tropicalis (424 read sets, MYA-3404 reference). Builds manuscript Figure 6
# and Tables 6-8. Reads existing outputs only; no pipeline is re-run here.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel); library(scales)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
sd   <- "/tmp/claude-1000/-mnt-HuangLab-26T-HuangLab-Data-Tar-Nine-Ta/6be50fd8-1892-479e-8dc9-2e9aa7ec1870/scratchpad"
td   <- file.path(base, "doc/Aspergillus/target_discovery_2026-08/tables")
pngdir <- file.path(base, "doc/manuscript_2026-10/figures"); figdir <- file.path(base, "output/figure")
theme_ms <- theme_classic(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 11), axis.text = element_text(colour = "black"),
        legend.key.size = unit(3, "mm"), strip.background = element_blank())
col_good <- "#1b7a5a"; col_bad <- "#b2182b"
save_fig <- function(p, n, w, h) {
  ggsave(file.path(figdir, paste0(n, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(pngdir, paste0(n, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}

af <- fread(file.path(td, "constraint_all_genes.csv"))
ct <- fread(file.path(base, "analysis/candida/constraint/ctrop_constraint_per_gene.tsv"))
om <- fread(file.path(base, "analysis/crossgenus/ortho_map.tsv"))

af_ok <- af[mean_cov_genomes >= 250 & !is.na(piNpiS_shrunk) & (syn_seg + nonsyn_seg) >= 3]
# Primary estimator for C. tropicalis is the RAW ratio on lineage representatives:
# this collection carries a median of ~108 segregating sites per gene, so the raw
# ratio is stable, whereas the genome-mean shrinkage used for the sparser
# Aspergillus data compresses these estimates and erases the ribosomal-protein
# signal (see analysis_code/ctrop_constraint.py).
ct_ok <- ct[n_iso_callable >= 300 & !is.na(piNpiS_dedup) & is.finite(piNpiS_dedup) & (syn_seg + nonsyn_seg) >= 3]
cat(sprintf("evaluable: Af %d/%d; Ct %d/%d genes\n", nrow(af_ok), nrow(af), nrow(ct_ok), nrow(ct)))
cat(sprintf("median piN/piS: Af %.3f (shrunk); Ct raw %.3f, raw-dedup %.3f, shrunk %.3f\n",
            median(af_ok$piNpiS_shrunk), median(ct_ok$piNpiS, na.rm = TRUE),
            median(ct_ok$piNpiS_dedup), median(ct_ok$piNpiS_shrunk)))
cat(sprintf("ribosomal proteins: Af %.3f vs %.3f genome; Ct %.3f vs %.3f genome\n",
            median(af_ok[grepl("ribosomal protein", product, ignore.case = TRUE), piNpiS_shrunk]),
            median(af_ok$piNpiS_shrunk),
            median(ct_ok[grepl("ribosomal protein", product, ignore.case = TRUE), piNpiS_dedup]),
            median(ct_ok$piNpiS_dedup)))
cat(sprintf("median segregating sites/gene: Af %.0f; Ct %.0f\n",
            median(af_ok$syn_seg + af_ok$nonsyn_seg), median(ct_ok$syn_seg + ct_ok$nonsyn_seg)))

m <- merge(om[ct_locus != "", .(af_protein, af_locus, ct_protein, ct_locus, ct_pident, calb_acc, calb_pident)],
           af_ok[, .(af_protein = protein_id, af_pi = piNpiS_shrunk, product, has_human_homolog)], by = "af_protein")
m <- merge(m, ct_ok[, .(ct_protein = protein_id, ct_pi = piNpiS_dedup, ct_pi_shrunk = piNpiS_shrunk_dedup, ct_product = product)], by = "ct_protein")
rho <- cor(m$af_pi, m$ct_pi, method = "spearman")
ctest <- cor.test(m$af_pi, m$ct_pi, method = "spearman")
cat(sprintf("n pairs = %d; Spearman rho = %.3f (p = %.2g); with shrunk Ct estimator %.3f\n",
            nrow(m), rho, ctest$p.value, cor(m$af_pi, m$ct_pi_shrunk, method = "spearman")))

m[, af_pct := 100 * (1 - frank(af_pi) / .N)]
m[, ct_pct := 100 * (1 - frank(ct_pi) / .N)]
q_af <- quantile(m$af_pi, .25); q_ct <- quantile(m$ct_pi, .25)
m[, quad := fifelse(af_pi <= q_af & ct_pi <= q_ct, "constrained in both",
             fifelse(af_pi <= q_af, "A. fumigatus only",
             fifelse(ct_pi <= q_ct, "C. tropicalis only", "neither")))]
print(m[, .N, by = quad][order(-N)])

# ---- ribosomal control: should be constrained in both ------------------------
m[, ribo := grepl("ribosomal protein", product, ignore.case = TRUE)]
cat(sprintf("ribosomal proteins (n=%d): constrained in both %.0f%% vs %.0f%% of the rest\n",
            m[ribo == TRUE, .N], 100 * m[ribo == TRUE, mean(quad == "constrained in both")],
            100 * m[ribo == FALSE, mean(quad == "constrained in both")]))

# ---- landmarks, including orthologue groups that are not one-to-one ----------
marks <- data.table(
  af_locus = c("AFUA_6G12400","AFUA_5G04170","AFUA_8G05210","AFUA_6G03540","AFUA_3G10310",
               "AFUA_4G13510","AFUA_2G05060","AFUA_1G13300","AFUA_7G03740"),
  name = c("fks1 / FKS1","hsp90 / HSP90","SAM-MTase (rank 1)","malate synthase","acetolactate synthase",
           "isocitrate lyase","alternative oxidase","GTP cyclohydrolase II / RIB1","cyp51B / ERG11"))
mk <- merge(marks, m, by = "af_locus")
print(mk[, .(name, af_locus, ct_locus, ident = ct_pident, af_pi = round(af_pi,3), af_pct = round(af_pct,1),
             ct_pi = round(ct_pi,3), ct_pct = round(ct_pct,1))])

f6a <- ggplot(m, aes(af_pct, ct_pct)) +
  geom_hline(yintercept = 75, linetype = "dotted", colour = "grey55") +
  geom_vline(xintercept = 75, linetype = "dotted", colour = "grey55") +
  geom_point(size = 0.3, alpha = 0.18, colour = "grey35") +
  geom_point(data = mk, colour = col_bad, size = 1.8) +
  geom_text_repel(data = mk, aes(label = name), size = 2.3, colour = col_bad,
                  min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf, seed = 1) +
  annotate("text", x = 2, y = 3, hjust = 0, size = 2.5,
           label = sprintf("rho = %.2f, n = %s", rho, comma(nrow(m)))) +
  scale_x_continuous(breaks = seq(0, 100, 25)) + scale_y_continuous(breaks = seq(0, 100, 25)) +
  labs(x = expression(italic("A. fumigatus")~"constraint percentile"),
       y = expression(italic("C. tropicalis")~"constraint percentile")) + theme_ms

# ---- candidates versus background -------------------------------------------
cand <- fread(file.path(td, "candidate_targets_ranked.csv"))
m[, is_cand := af_locus %in% cand$locus_tag]
orth_any <- om[ct_locus != "", unique(af_locus)]
cat(sprintf("candidates with a 1:1 C. tropicalis orthologue: %d/%d (%.0f%%); genome baseline %.0f%%\n",
            sum(cand$locus_tag %in% orth_any), nrow(cand),
            100 * mean(cand$locus_tag %in% orth_any), 100 * mean(af$locus_tag %in% orth_any)))
w <- wilcox.test(ct_pct ~ is_cand, data = m)
cat(sprintf("Ct constraint percentile, candidates %.1f vs others %.1f (median), Wilcoxon p = %.3g\n",
            m[is_cand == TRUE, median(ct_pct)], m[is_cand == FALSE, median(ct_pct)], w$p.value))
afq <- m[af_pi <= q_af]
cat(sprintf("among genes already constrained in A. fumigatus (n=%d), %.0f%% are also in the C. tropicalis top quartile (chance = 25%%)\n",
            nrow(afq), 100 * mean(afq$ct_pi <= q_ct)))

grp_dt <- rbind(
  data.table(grp = "All\npairs", ct_pct = m$ct_pct),
  data.table(grp = "Af-constrained,\nhuman hom.", ct_pct = m[af_pi <= q_af & has_human_homolog == TRUE, ct_pct]),
  data.table(grp = "Af-constrained,\nno human hom.", ct_pct = m[af_pi <= q_af & has_human_homolog == FALSE, ct_pct]),
  data.table(grp = "Prioritised\ncandidates", ct_pct = m[is_cand == TRUE, ct_pct]))
grp_dt[, grp := factor(grp, levels = unique(grp))]
nlab <- grp_dt[, .(n = .N, med = median(ct_pct)), by = grp]
print(nlab)
wh <- wilcox.test(ct_pct ~ has_human_homolog, data = m[af_pi <= q_af])
cat(sprintf("among Af-constrained genes, Ct percentile with vs without human homologue: p = %.2g\n", wh$p.value))
f6b <- ggplot(grp_dt, aes(grp, ct_pct, fill = grp)) +
  geom_hline(yintercept = 50, linetype = "dotted", colour = "grey55") +
  geom_boxplot(outlier.size = 0.3, width = 0.6, colour = "grey30", linewidth = 0.3) +
  geom_text(data = nlab, aes(grp, 103, label = paste0("n=", n)), size = 2.1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("grey75", "#9ecae1", "#fdd49e", col_good), guide = "none") +
  labs(x = NULL, y = expression(italic("C. tropicalis")~"constraint percentile")) +
  theme_ms + theme(axis.text.x = element_text(size = 6.3))

# ---- resistance-side convergence at the azole target ------------------------
rm_ <- fread(file.path(base, "doc/Candida/tables/erg11_cyp51A_residue_map.csv"))
cta <- fread(file.path(base, "doc/Candida/tables/ctropicalis_analysis_table.csv"))[qc_pass == TRUE]
afg <- fread(file.path(base, "output/table/vrc_cyp51A_genotypes.csv"))
conv <- rbind(
  data.table(genus = "C. tropicalis (fluconazole)", allele = "ERG11 Y132F",
             carrier = cta$Y132F == 1, mic = cta$mic_used),
  data.table(genus = "C. tropicalis (fluconazole)", allele = "ERG11 G464S",
             carrier = cta$G464S == 1, mic = cta$mic_used),
  data.table(genus = "A. fumigatus (voriconazole)", allele = "cyp51A Y121F",
             carrier = grepl("Y121F", afg$aa_substitutions), mic = afg$mic_used),
  data.table(genus = "A. fumigatus (voriconazole)", allele = "cyp51A G448S",
             carrier = grepl("G448S", afg$aa_substitutions), mic = afg$mic_used))
conv[, grp := fifelse(carrier, "carrier", "non-carrier")]
lab <- conv[, .(n = .N, med = median(mic)), by = .(genus, allele, grp)]
print(lab)
conv[, allele := factor(allele, levels = c("ERG11 Y132F", "cyp51A Y121F", "ERG11 G464S", "cyp51A G448S"))]
f6c <- ggplot(conv, aes(allele, mic, colour = grp, fill = grp)) +
  geom_boxplot(outlier.shape = NA, width = .6, position = position_dodge(.75), alpha = .15, linewidth = .3) +
  geom_point(position = position_jitterdodge(jitter.width = .14, dodge.width = .75), size = .7, alpha = .6) +
  scale_fill_manual(values = c(carrier = col_bad, `non-carrier` = "grey50"), guide = "none") +
  scale_y_continuous(trans = "log2", breaks = 2^seq(-2, 10, 2)) +
  scale_colour_manual(values = c(carrier = col_bad, `non-carrier` = "grey50"), name = NULL) +
  facet_wrap(~ genus, scales = "free", nrow = 1) +
  labs(x = NULL, y = "Azole MIC (mg/L)") + theme_ms +
  theme(axis.text.x = element_text(size = 7), legend.position = "bottom")

f6 <- (f6a | f6b) / f6c + plot_layout(heights = c(1.15, 1)) + plot_annotation(tag_levels = "A")
save_fig(f6, "ms_fig6_crossgenus", 7.2, 6.6)

# ---- tables -----------------------------------------------------------------
fwrite(mk[, .(Gene = name, Af_locus = af_locus, Ct_locus = ct_locus, identity_pct = ct_pident,
              Af_piNpiS = round(af_pi, 3), Af_pctile = round(af_pct, 1),
              Ct_piNpiS = round(ct_pi, 3), Ct_pctile = round(ct_pct, 1))],
       file.path(base, "output/table/ms_table6_crossgenus_landmarks.csv"))
t7 <- m[is_cand == TRUE][order(ct_pi)][, .(Af_locus = af_locus, Product = product,
        Ct_orthologue = ct_locus, identity_pct = ct_pident, Af_piNpiS = round(af_pi, 3),
        Af_pctile = round(af_pct, 1), Ct_piNpiS = round(ct_pi, 3), Ct_pctile = round(ct_pct, 1), Quadrant = quad)]
fwrite(t7, file.path(base, "output/table/ms_table7_crossgenus_candidates.csv"))
fwrite(rm_, file.path(base, "output/table/ms_table8_residue_convergence.csv"))
cat("\ncandidates constrained in both genera:\n"); print(head(t7[Quadrant == "constrained in both"], 20))
saveRDS(m, file.path(sd, "crossgenus_merged.rds"))

asmdir <- file.path(base, "doc/manuscript_2026-10/figures_mSystems"); dir.create(asmdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(asmdir, "FIG3.tiff"), f6, width = 7, height = 6.6, dpi = 300, compression = "lzw", bg = "white")
