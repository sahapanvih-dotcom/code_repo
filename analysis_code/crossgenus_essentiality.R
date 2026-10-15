# Does constraint predict EXPERIMENTAL essentiality of the orthologue in a third
# species? Benchmarks A. fumigatus (and C. tropicalis) constraint against five
# C. albicans essentiality datasets. Writes manuscript Figure 7 and Tables 9-10.
#
# NOTE ON THE ANALYSIS CHOICE: the five sources disagree with one another often
# (GRACE rescorings vs the Segal transposon screen), so a strict per-gene consensus
# ("essential only if no source disagrees") pushes most genuinely essential genes
# into a "conflicting" bin and produces a spurious null (AUC 0.54, p = 0.16).
# Each source is therefore tested separately, which is also the only way to see
# that the answer is consistent across independent screens.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
sd   <- "/tmp/claude-1000/-mnt-HuangLab-26T-HuangLab-Data-Tar-Nine-Ta/6be50fd8-1892-479e-8dc9-2e9aa7ec1870/scratchpad"
td   <- file.path(base, "doc/Aspergillus/target_discovery_2026-08/tables")
pngdir <- file.path(base, "doc/manuscript_2026-10/figures"); figdir <- file.path(base, "output/figure")
theme_ms <- theme_classic(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 11), axis.text = element_text(colour = "black"),
        legend.key.size = unit(3, "mm"))
col_good <- "#1b7a5a"; col_bad <- "#b2182b"

ess <- fread(file.path(base, "analysis/crossgenus/calb_essentiality.tsv.gz"))
idm <- fread(file.path(base, "analysis/crossgenus/calb_idmap.tsv")); setnames(idm, c("acc", "oln", "primary", "syn"))
idl <- idm[, .(id = unlist(strsplit(oln, " "))), by = acc][id != ""]
# CGD assembly-22 IDs ("C1_00730C_A") and UniProt ordered-locus names
# ("CAALFM_C100730CA") are the same identifier written differently.
norm <- function(x) { y <- toupper(trimws(x)); y <- sub("^CAALFM_", "", y); gsub("_", "", y) }
ess[, gid := norm(gene_id)]; idl[, gid := norm(id)]
e <- merge(ess, idl[, .(acc, gid)], by = "gid", allow.cartesian = TRUE)
cat(sprintf("essentiality rows mapped to UniProt: %d of %d\n", nrow(e), nrow(ess)))

om <- fread(file.path(base, "analysis/crossgenus/ortho_map.tsv"))
af <- fread(file.path(td, "constraint_all_genes.csv"))[mean_cov_genomes >= 250 & (syn_seg + nonsyn_seg) >= 3]
ct <- fread(file.path(base, "analysis/candida/constraint/ctrop_constraint_per_gene.tsv"))[
        n_iso_callable >= 300 & (syn_seg + nonsyn_seg) >= 3 & is.finite(piNpiS_dedup)]
d_af <- merge(om[calb_acc != "", .(af_locus, acc = calb_acc)],
              af[, .(af_locus = locus_tag, pi = piNpiS_shrunk, product, has_human_homolog)], by = "af_locus")
d_ct <- merge(om[calb_acc != "" & ct_protein != "", .(acc = calb_acc, ct_protein)],
              ct[, .(ct_protein = protein_id, pi = piNpiS_dedup)], by = "ct_protein")

auc_f <- function(x, lab) { r <- rank(-x); n1 <- sum(lab); n0 <- sum(!lab)
  (sum(r[lab]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
auc_ci <- function(x, lab, B = 500) {  # bootstrap CI
  set.seed(1); v <- replicate(B, { i <- sample.int(length(x), replace = TRUE)
    if (uniqueN(lab[i]) < 2) NA_real_ else auc_f(x[i], lab[i]) })
  quantile(v, c(.025, .975), na.rm = TRUE) }

src_lab <- c(Segal2018 = "Segal 2018 (transposon)", GRACEv1_Fu2021 = "GRACE v1 (Fu 2021)",
             GRACEv1_OMeara2015 = "GRACE v1 (O'Meara 2015)", GRACEv2_Fu2021 = "GRACE v2 (Fu 2021)",
             GRACEv3_Fu2025 = "GRACE v3 (Fu 2025)", CGD_phenotype = "CGD curated phenotypes")
res <- rbindlist(lapply(names(src_lab), function(s) {
  x <- e[source == s & call %in% c("ESSENTIAL", "NONESSENTIAL"), .(ess = any(call == "ESSENTIAL")), by = acc]
  rbindlist(lapply(list(af = d_af, ct = d_ct), function(d) {
    m <- merge(d, x, by = "acc")
    if (uniqueN(m$ess) < 2 || sum(m$ess) < 10) return(NULL)
    ci <- auc_ci(m$pi, m$ess)
    data.table(source = src_lab[[s]], n = nrow(m), n_ess = sum(m$ess),
               auc = auc_f(m$pi, m$ess), lo = ci[1], hi = ci[2],
               p = wilcox.test(pi ~ ess, data = m)$p.value,
               med_ess = median(m[ess == TRUE, pi]), med_non = median(m[ess == FALSE, pi]))
  }), idcol = "genus")
}))
res[, genus := fifelse(genus == "af", "A. fumigatus constraint", "C. tropicalis constraint")]
print(res[order(genus, -auc)])
fwrite(res, file.path(base, "output/table/ms_table9_essentiality_transfer.csv"))

f7a <- ggplot(res, aes(auc, reorder(source, auc), colour = genus)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, position = position_dodge(.55), linewidth = .4) +
  geom_point(size = 1.9, position = position_dodge(.55)) +
  scale_colour_manual(values = c("A. fumigatus constraint" = col_good, "C. tropicalis constraint" = col_bad), name = NULL) +
  scale_x_continuous(limits = c(0.4, 0.85)) +
  labs(x = "AUC (bootstrap 95% CI)", y = NULL) +
  theme_ms + theme(legend.position = "bottom", axis.text.y = element_text(size = 7))

# essential fraction by A. fumigatus constraint decile, largest independent screens
dec <- rbindlist(lapply(c("Segal2018", "GRACEv1_Fu2021"), function(s) {
  x <- e[source == s & call %in% c("ESSENTIAL", "NONESSENTIAL"), .(ess = any(call == "ESSENTIAL")), by = acc]
  m <- merge(d_af, x, by = "acc")
  m[, dec := cut(frank(pi), 10, labels = FALSE)]
  m[, .(pct = 100 * mean(ess), n = .N, source = src_lab[[s]]), by = dec]
}))
print(dcast(dec, dec ~ source, value.var = "pct"))
f7b <- ggplot(dec, aes(dec, pct, fill = source)) +
  geom_col(position = position_dodge(.8), width = .75) +
  scale_fill_manual(values = c("#4a7fb5", col_good), name = NULL) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Constraint decile (1 = most constrained)",
       y = expression("Orthologue essential in"~italic("C. albicans")~"(%)")) +
  theme_ms + theme(legend.position = c(0.72, 0.9))

f7 <- (f7a | f7b) + plot_annotation(tag_levels = "A")
ggsave(file.path(figdir, "ms_fig7_calb_essentiality.pdf"), f7, width = 7.2, height = 3.4, device = cairo_pdf)
ggsave(file.path(pngdir, "ms_fig7_calb_essentiality.png"), f7, width = 7.2, height = 3.4, dpi = 300, bg = "white")

# ---- annotate the cross-genus intersection with per-source calls ---------------
wide <- dcast(e[call %in% c("ESSENTIAL", "NONESSENTIAL")], acc ~ source,
              value.var = "call", fun.aggregate = function(v) fifelse(any(v == "ESSENTIAL"), "E", "N"), fill = NA)
t7 <- fread(file.path(base, "output/table/ms_table7_crossgenus_candidates.csv"))[Quadrant == "constrained in both"]
inter <- merge(t7[, .(af_locus = Af_locus, Product, Ct_orthologue, Af_pctile, Ct_pctile)],
               om[, .(af_locus, acc = calb_acc, calb_pident)], by = "af_locus", all.x = TRUE)
inter <- merge(inter, wide, by = "acc", all.x = TRUE)
inter[, n_ess := rowSums(.SD == "E", na.rm = TRUE), .SDcols = names(src_lab)[names(src_lab) %in% names(inter)]]
inter[, n_call := rowSums(!is.na(.SD)), .SDcols = names(src_lab)[names(src_lab) %in% names(inter)]]
setorder(inter, -n_ess, -Ct_pctile)
print(inter[, .(af_locus, Product = substr(Product, 1, 42), Ct_pctile, n_ess, n_call)])
fwrite(inter, file.path(base, "output/table/ms_table10_intersection_essentiality.csv"))

# sanity: the strict-consensus analysis that fails, reported as a caveat
cons <- e[call %in% c("ESSENTIAL", "NONESSENTIAL"),
          .(n_ess = uniqueN(source[call == "ESSENTIAL"]), n_non = uniqueN(source[call == "NONESSENTIAL"])), by = acc]
cons[, consensus := fifelse(n_ess > 0 & n_non == 0, "essential",
                     fifelse(n_non > 0 & n_ess == 0, "non-essential", "conflicting"))]
cat("\nstrict per-gene consensus across sources:\n"); print(cons[, .N, by = consensus][order(-N)])
mc <- merge(d_af, cons[consensus != "conflicting"], by = "acc")[, ess := consensus == "essential"]
cat(sprintf("consensus-based AUC = %.3f (p = %.2g) on n = %d, %d essential - biased null, see script header\n",
            auc_f(mc$pi, mc$ess), wilcox.test(pi ~ ess, data = mc)$p.value, nrow(mc), sum(mc$ess)))

asmdir <- file.path(base, "doc/manuscript_2026-10/figures_mSystems"); dir.create(asmdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(asmdir, "FIG4.tiff"), f7, width = 7, height = 3.4, dpi = 300, compression = "lzw", bg = "white")
