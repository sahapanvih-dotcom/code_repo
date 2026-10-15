#!/usr/bin/env Rscript
# rev_b4_02a_coverage.R
# Per-sequence alignment coverage for the PF00925 (RIB1) family alignment, plus the
# Pfam domain architecture table. R is used for all summarising/statistics per lab
# convention; ape reads the alignment.
suppressPackageStartupMessages(library(ape))

BASE <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
OG   <- file.path(BASE, "analysis/orthogroup")

tips <- read.delim(file.path(OG, "fasta/rev_b4_tiplabels.tsv"), stringsAsFactors = FALSE)

aln_files <- c(
  mafft_linsi_untrimmed = file.path(OG, "aln/rev_b4_family_mafft_linsi.faa"),
  trim_automated1       = file.path(OG, "trim/rev_b4_mafft_automated1.faa"),
  trim_gt08             = file.path(OG, "trim/rev_b4_mafft_gt08.faa"),
  trim_gt03             = file.path(OG, "trim/rev_b4_mafft_gt03.faa")
)
aln_files <- aln_files[file.exists(aln_files)]

cov_list <- lapply(names(aln_files), function(nm) {
  a <- read.FASTA(aln_files[[nm]], type = "AA")
  m <- as.character(as.matrix(a))
  ncol_aln <- ncol(m)
  gapchars <- c("-", "?", ".", "x")           # ape lower-cases AA matrices
  occ <- apply(m, 1, function(r) sum(!(tolower(r) %in% gapchars)))
  data.frame(alignment = nm, tip = rownames(m), aln_ncol = ncol_aln,
             residues_in_aln = as.integer(occ),
             aln_coverage_frac = round(occ / ncol_aln, 4),
             stringsAsFactors = FALSE)
})
cov <- do.call(rbind, cov_list)

# ---- Pfam domain architecture from the step-1 hmmscan domtblout
dom <- read.table(file.path(OG, "pfam/rev_b4_candidates_pfam.domtbl"),
                  comment.char = "#", fill = TRUE, stringsAsFactors = FALSE)
dom <- dom[, c(1, 2, 4, 6, 7, 12, 18, 19)]
names(dom) <- c("pfam_name", "pfam_acc", "acc", "aa_len", "seq_evalue",
                "dom_ievalue", "ali_from", "ali_to")
dom <- dom[dom$pfam_name %in% c("GTP_cyclohydro2", "GTP_CH_N", "DHDPS"), ]
write.csv(dom[order(dom$acc, dom$ali_from), ],
          file.path(OG, "pfam/rev_b4_domain_coordinates.csv"), row.names = FALSE)

arch <- tapply(dom$pfam_name, dom$acc, function(x) paste(sort(unique(x)), collapse = "+"))
# best (smallest) PF00925 domain E-value per sequence: how intact the catalytic domain is
pf925 <- dom[dom$pfam_name == "GTP_cyclohydro2", ]
best925 <- tapply(pf925$seq_evalue, pf925$acc, min)

tips$domain_architecture <- arch[tips$acc]
tips$PF00925_seq_evalue  <- signif(as.numeric(best925[tips$acc]), 3)
tips$has_PF12471_GTP_CH_N <- grepl("GTP_CH_N", tips$domain_architecture)

out <- merge(tips, cov[cov$alignment == names(aln_files)[1],
                       c("tip", "aln_ncol", "residues_in_aln", "aln_coverage_frac")],
             by = "tip", all.x = TRUE)
names(out)[names(out) == "aln_ncol"]          <- "untrimmed_aln_ncol"
names(out)[names(out) == "residues_in_aln"]   <- "untrimmed_residues"
names(out)[names(out) == "aln_coverage_frac"] <- "untrimmed_coverage"

if ("trim_automated1" %in% cov$alignment) {
  t1 <- cov[cov$alignment == "trim_automated1", c("tip", "residues_in_aln", "aln_coverage_frac", "aln_ncol")]
  names(t1) <- c("tip", "trimmed_residues", "trimmed_coverage", "trimmed_aln_ncol")
  out <- merge(out, t1, by = "tip", all.x = TRUE)
}

out <- out[order(out$has_PF12471_GTP_CH_N, out$species), ]
write.csv(out, file.path(OG, "rev_b4_family_coverage_domains.csv"), row.names = FALSE)
write.csv(cov, file.path(OG, "rev_b4_alignment_coverage_all.csv"), row.names = FALSE)

cat("alignment columns:\n"); print(unique(cov[, c("alignment", "aln_ncol")]))
cat("\nmean/min coverage (untrimmed L-INS-i):",
    round(mean(out$untrimmed_coverage), 3), "/", round(min(out$untrimmed_coverage), 3), "\n")
cat("\ndomain architecture counts:\n"); print(table(out$domain_architecture))
