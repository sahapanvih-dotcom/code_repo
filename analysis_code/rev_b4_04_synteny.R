#!/usr/bin/env Rscript
# rev_b4_04_synteny.R
# Revision analysis B4, step 4: conserved-synteny evidence for the RIB1 orthology call.
#
# For each of the four loci of interest the +/-5 flanking protein-coding genes are read
# straight out of the RefSeq GFF annotations, their proteins are extracted, and the
# A. fumigatus neighbourhoods are compared with the Candida / S. cerevisiae RIB1
# neighbourhoods by reciprocal blastp. "Shared neighbours" = flanking genes whose best
# blastp hit in the other neighbourhood is reciprocal (E <= 1e-5).
#
# R does the GFF parsing, the window logic and all summarising; blastp (2.12.0+) and
# seqkit do the sequence steps, which is what they are for.
suppressPackageStartupMessages(library(ape))

BASE <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
OG   <- file.path(BASE, "analysis/orthogroup")
SYN  <- file.path(OG, "synteny")
TAB  <- file.path(BASE, "output/table")
dir.create(SYN, showWarnings = FALSE, recursive = TRUE)
FLANK <- 5

loci <- data.frame(
  label   = c("Af_AFUA_1G13300", "Af_AFUA_2G01220", "Ct_CTRG_00799_RIB1",
              "Ca_RIB1", "Sc_RIB1"),
  species = c("Aspergillus_fumigatus_Af293", "Aspergillus_fumigatus_Af293",
              "Candida_tropicalis_MYA3404", "Candida_albicans_SC5314",
              "Saccharomyces_cerevisiae_S288C"),
  protein = c("XP_077660111.1", "XP_749258.2", "XP_002546018.1",
              "XP_715792.1", "NP_009520.1"),
  gff = c(file.path(BASE, "raw_data/coverage_pass/Aspergillus_fumigatus/GCF_000002655.1/genomic.gff"),
          file.path(BASE, "raw_data/coverage_pass/Aspergillus_fumigatus/GCF_000002655.1/genomic.gff"),
          file.path(BASE, "analysis/candida/ref/genomic.gff"),
          file.path(OG, "ref/calb_sc5314.gff"),
          file.path(OG, "ref/scer_s288c.gff")),
  faa = c(file.path(BASE, "raw_data/coverage_pass/Aspergillus_fumigatus/GCF_000002655.1/protein.faa"),
          file.path(BASE, "raw_data/coverage_pass/Aspergillus_fumigatus/GCF_000002655.1/protein.faa"),
          file.path(BASE, "analysis/candida/ref/protein.faa"),
          file.path(OG, "ref/calb_sc5314.faa"),
          file.path(OG, "ref/scer_s288c.faa")),
  stringsAsFactors = FALSE)
stopifnot(all(file.exists(loci$gff)), all(file.exists(loci$faa)))

attr_get <- function(x, key) {
  pat <- paste0("(^|;)", key, "=")
  pos <- regexpr(pat, x)
  out <- rep(NA_character_, length(x))
  ok  <- pos > 0
  if (any(ok)) {
    tail_ <- substring(x[ok], pos[ok] + attr(pos, "match.length")[ok])
    out[ok] <- sub(";.*$", "", tail_)
  }
  out
}

read_cds <- function(gff) {
  g <- read.delim(gff, header = FALSE, comment.char = "#", quote = "",
                  col.names = c("seqid","src","type","start","end","score","strand","phase","attr"),
                  stringsAsFactors = FALSE)
  g <- g[g$type == "CDS", ]
  g$protein_id <- attr_get(g$attr, "protein_id")
  g$locus_tag  <- attr_get(g$attr, "locus_tag")
  g$gene       <- attr_get(g$attr, "gene")
  g$product    <- attr_get(g$attr, "product")
  g <- g[!is.na(g$protein_id), ]
  g <- g[!duplicated(g$protein_id), ]
  g[order(g$seqid, g$start), c("seqid","start","end","strand","protein_id","locus_tag","gene","product")]
}

windows <- list()
for (i in seq_len(nrow(loci))) {
  cds <- read_cds(loci$gff[i])
  hit <- which(cds$protein_id == loci$protein[i])
  if (!length(hit)) { warning("anchor not found: ", loci$label[i]); next }
  hit <- hit[1]
  same <- which(cds$seqid == cds$seqid[hit])
  pos  <- match(hit, same)
  idx  <- same[seq(max(1, pos - FLANK), min(length(same), pos + FLANK))]
  w <- cds[idx, ]
  w$anchor_label  <- loci$label[i]
  w$anchor_species <- loci$species[i]
  w$offset <- seq_along(idx) - match(hit, idx)
  w$is_anchor <- w$offset == 0
  windows[[loci$label[i]]] <- w
}
allw <- do.call(rbind, windows)
write.csv(allw, file.path(SYN, "rev_b4_synteny_windows.csv"), row.names = FALSE)

# extract each window's proteins
for (nm in names(windows)) {
  i <- match(nm, loci$label)
  ids <- file.path(SYN, paste0(nm, ".ids"))
  writeLines(windows[[nm]]$protein_id, ids)
  system2("seqkit", c("grep", "-f", shQuote(ids), shQuote(loci$faa[i])),
          stdout = file.path(SYN, paste0(nm, ".faa")), stderr = file.path(SYN, paste0(nm, ".seqkit.log")))
}

# reciprocal blastp between every A. fumigatus window and every yeast window
bl <- function(q, s, out) {
  system2("blastp", c("-query", shQuote(q), "-subject", shQuote(s), "-evalue", "1e-5",
                      "-outfmt", shQuote("6 qseqid sseqid pident length qlen evalue bitscore"),
                      "-max_hsps", "1"), stdout = out, stderr = paste0(out, ".log"))
  if (file.info(out)$size == 0) return(NULL)
  d <- read.delim(out, header = FALSE,
                  col.names = c("q","s","pident","len","qlen","evalue","bits"))
  d[order(d$q, -d$bits), ]
}
best <- function(d) if (is.null(d)) NULL else d[!duplicated(d$q), ]

pairs <- expand.grid(af = c("Af_AFUA_1G13300", "Af_AFUA_2G01220"),
                     yl = c("Ct_CTRG_00799_RIB1", "Ca_RIB1", "Sc_RIB1"),
                     stringsAsFactors = FALSE)
syn <- lapply(seq_len(nrow(pairs)), function(k) {
  a <- pairs$af[k]; b <- pairs$yl[k]
  fa <- file.path(SYN, paste0(a, ".faa")); fb <- file.path(SYN, paste0(b, ".faa"))
  f <- best(bl(fa, fb, file.path(SYN, paste0(a, "__vs__", b, ".tsv"))))
  r <- best(bl(fb, fa, file.path(SYN, paste0(b, "__vs__", a, ".tsv"))))
  wa <- windows[[a]]; wb <- windows[[b]]
  anchor_a <- wa$protein_id[wa$is_anchor]; anchor_b <- wb$protein_id[wb$is_anchor]
  recip <- character(0); det <- NULL
  if (!is.null(f) && !is.null(r)) {
    for (j in seq_len(nrow(f))) {
      q <- f$q[j]; s <- f$s[j]
      back <- r$s[r$q == s]
      if (length(back) && back[1] == q) {
        recip <- c(recip, q)
        det <- rbind(det, data.frame(af_window = a, yeast_window = b,
                                     af_protein = q, af_offset = wa$offset[match(q, wa$protein_id)],
                                     af_locus = wa$locus_tag[match(q, wa$protein_id)],
                                     af_product = wa$product[match(q, wa$protein_id)],
                                     yeast_protein = s, yeast_offset = wb$offset[match(s, wb$protein_id)],
                                     yeast_gene = wb$gene[match(s, wb$protein_id)],
                                     yeast_product = wb$product[match(s, wb$protein_id)],
                                     pident = f$pident[j], evalue = f$evalue[j], bits = f$bits[j],
                                     stringsAsFactors = FALSE))
      }
    }
  }
  anchor_recip <- anchor_a %in% recip &&
    !is.null(det) && any(det$af_protein == anchor_a & det$yeast_protein == anchor_b)
  flank_recip <- setdiff(recip, anchor_a)
  list(summary = data.frame(af_window = a, yeast_window = b,
                            n_genes_af_window = nrow(wa), n_genes_yeast_window = nrow(wb),
                            anchor_is_reciprocal_best = anchor_recip,
                            n_reciprocal_pairs_total = length(recip),
                            n_reciprocal_flanking_pairs = length(flank_recip),
                            stringsAsFactors = FALSE),
       detail = det)
})
syn_summary <- do.call(rbind, lapply(syn, `[[`, "summary"))
syn_detail  <- do.call(rbind, lapply(syn, `[[`, "detail"))

write.csv(syn_summary, file.path(TAB, "rev_b4_rib1_synteny_summary.csv"), row.names = FALSE)
if (!is.null(syn_detail))
  write.csv(syn_detail, file.path(TAB, "rev_b4_rib1_synteny_shared_neighbours.csv"), row.names = FALSE)

cat("\n== synteny windows (+/-", FLANK, "genes) ==\n")
print(as.data.frame(table(allw$anchor_label)))
cat("\n== neighbourhood comparison ==\n"); print(syn_summary)
cat("\n== reciprocal pairs ==\n")
if (!is.null(syn_detail)) print(syn_detail[, c("af_window","yeast_window","af_locus","af_offset",
                                               "yeast_gene","yeast_offset","pident","bits")]) else cat("none\n")
cat("\n== anchor neighbourhoods ==\n")
print(allw[, c("anchor_label","offset","seqid","start","end","strand","locus_tag","gene","product")],
      row.names = FALSE)

# ======================================================================================
# Part 2: synteny at phylogenetic distances where microsynteny is actually informative.
# Aspergillus vs Candida/Saccharomyces microsynteny is NOT expected to survive (reported
# above as 0 shared flanking neighbours), so the informative comparisons are
#   (i)  Candida vs Candida vs Saccharomyces around RIB1  - does CTRG_00799 sit in the
#        same neighbourhood as the C. albicans RIB1 that was validated?
#   (ii) Aspergillus vs Aspergillus around each class     - are the AFUA_1G13300-type and
#        AFUA_2G01220-type loci two distinct, independently conserved genomic positions?
# ======================================================================================
loci2 <- data.frame(
  label   = c("Alen_classI_XP033411727", "Alen_classII_XP033412737",
              "Anid_classI_XP658274",    "Anid_classII_XP050468011",
              "Anig_classI_XP059604098", "Anig_classII_XP001402428"),
  species = c("Aspergillus_lentulus","Aspergillus_lentulus",
              "Aspergillus_nidulans","Aspergillus_nidulans",
              "Aspergillus_niger","Aspergillus_niger"),
  protein = c("XP_033411727.1","XP_033412737.1","XP_658274.1",
              "XP_050468011.1","XP_059604098.1","XP_001402428.1"),
  acc     = c("GCF_010724455.1","GCF_010724455.1","GCF_000011425.1",
              "GCF_000011425.1","GCF_000002855.4","GCF_000002855.4"),
  sp_dir  = c("Aspergillus_lentulus","Aspergillus_lentulus","Aspergillus_nidulans",
              "Aspergillus_nidulans","Aspergillus_niger","Aspergillus_niger"),
  stringsAsFactors = FALSE)
loci2$gff <- file.path(BASE, "raw_data/coverage_pass", loci2$sp_dir, loci2$acc, "genomic.gff")
loci2$faa <- file.path(BASE, "raw_data/coverage_pass", loci2$sp_dir, loci2$acc, "protein.faa")
loci2 <- loci2[file.exists(loci2$gff) & file.exists(loci2$faa), ]

for (i in seq_len(nrow(loci2))) {
  cds <- read_cds(loci2$gff[i]); hit <- which(cds$protein_id == loci2$protein[i])
  if (!length(hit)) { warning("anchor not found: ", loci2$label[i]); next }
  hit <- hit[1]; same <- which(cds$seqid == cds$seqid[hit]); pos <- match(hit, same)
  idx <- same[seq(max(1, pos - FLANK), min(length(same), pos + FLANK))]
  w <- cds[idx, ]; w$anchor_label <- loci2$label[i]; w$anchor_species <- loci2$species[i]
  w$offset <- seq_along(idx) - match(hit, idx); w$is_anchor <- w$offset == 0
  windows[[loci2$label[i]]] <- w
  ids <- file.path(SYN, paste0(loci2$label[i], ".ids"))
  writeLines(w$protein_id, ids)
  system2("seqkit", c("grep", "-f", shQuote(ids), shQuote(loci2$faa[i])),
          stdout = file.path(SYN, paste0(loci2$label[i], ".faa")),
          stderr = file.path(SYN, paste0(loci2$label[i], ".seqkit.log")))
}

compare_windows <- function(a, b) {
  fa <- file.path(SYN, paste0(a, ".faa")); fb <- file.path(SYN, paste0(b, ".faa"))
  if (!file.exists(fa) || !file.exists(fb)) return(NULL)
  f <- best(bl(fa, fb, file.path(SYN, paste0(a, "__vs__", b, ".tsv"))))
  r <- best(bl(fb, fa, file.path(SYN, paste0(b, "__vs__", a, ".tsv"))))
  wa <- windows[[a]]; wb <- windows[[b]]
  anchor_a <- wa$protein_id[wa$is_anchor]; anchor_b <- wb$protein_id[wb$is_anchor]
  recip <- character(0); anchor_ok <- FALSE
  if (!is.null(f) && !is.null(r)) for (j in seq_len(nrow(f))) {
    q <- f$q[j]; sj <- f$s[j]; back <- r$s[r$q == sj]
    if (length(back) && back[1] == q) {
      recip <- c(recip, q)
      if (q == anchor_a && sj == anchor_b) anchor_ok <- TRUE
    }
  }
  data.frame(window_A = a, window_B = b,
             anchor_is_reciprocal_best = anchor_ok,
             n_reciprocal_pairs_total = length(recip),
             n_reciprocal_flanking_pairs = length(setdiff(recip, anchor_a)),
             flanking_synteny_frac = round(length(setdiff(recip, anchor_a)) / (nrow(wa) - 1), 2),
             stringsAsFactors = FALSE)
}

cmp <- list(
  c("Ct_CTRG_00799_RIB1", "Ca_RIB1"),
  c("Ct_CTRG_00799_RIB1", "Sc_RIB1"),
  c("Ca_RIB1", "Sc_RIB1"),
  c("Af_AFUA_1G13300", "Alen_classI_XP033411727"),
  c("Af_AFUA_2G01220", "Alen_classII_XP033412737"),
  c("Af_AFUA_1G13300", "Anid_classI_XP658274"),
  c("Af_AFUA_2G01220", "Anid_classII_XP050468011"),
  c("Af_AFUA_1G13300", "Anig_classI_XP059604098"),
  c("Af_AFUA_2G01220", "Anig_classII_XP001402428"),
  c("Af_AFUA_1G13300", "Af_AFUA_2G01220"),
  c("Af_AFUA_1G13300", "Anid_classII_XP050468011"),
  c("Af_AFUA_2G01220", "Anid_classI_XP658274"))
syn2 <- do.call(rbind, lapply(cmp, function(z) compare_windows(z[1], z[2])))
syn2$comparison_type <- c("within Saccharomycotina (RIB1)", "within Saccharomycotina (RIB1)",
  "within Saccharomycotina (RIB1)", "Aspergillus class I vs class I",
  "Aspergillus class II vs class II", "Aspergillus class I vs class I",
  "Aspergillus class II vs class II", "Aspergillus class I vs class I",
  "Aspergillus class II vs class II", "class I vs class II (same genome)",
  "class I vs class II (cross-species, negative control)",
  "class II vs class I (cross-species, negative control)")[seq_len(nrow(syn2))]
write.csv(syn2, file.path(TAB, "rev_b4_rib1_synteny_within_lineage.csv"), row.names = FALSE)
allw2 <- do.call(rbind, windows)
write.csv(allw2, file.path(SYN, "rev_b4_synteny_windows.csv"), row.names = FALSE)
cat("\n== within-lineage synteny (the informative comparisons) ==\n"); print(syn2)

# ======================================================================================
# Part 3: the class-II locus is a conserved TWO-GENE CASSETTE, fused in some species.
# AFUA_2G01220 is immediately followed by AFUA_2G01230 (dihydrodipicolinate synthase
# family). Two class-II members found in step 1 are ~785 aa and carry PF00701 (DHDPS) in
# addition to PF12471+PF00925. If they are fusions of exactly that adjacent gene pair,
# the class-II locus is a Pezizomycotina-restricted gene module, which no Saccharomycotina
# genome has - independent of any tree.
# ======================================================================================
fus_ids <- c("XP_001216775.1", "XP_024676616.1")   # A. terreus, A. candidus class-II fusions
fus_faa <- file.path(SYN, "rev_b4_classII_fusions.faa")
writeLines(fus_ids, file.path(SYN, "rev_b4_classII_fusions.ids"))
fam_faa <- file.path(OG, "fasta/rev_b4_candidates.faa")
keep <- FALSE; out <- character(0)
for (ln in readLines(fam_faa)) {
  if (startsWith(ln, ">")) keep <- sub("^>([^ ]+).*$", "\\1", ln) %in% fus_ids
  if (keep) out <- c(out, ln)
}
writeLines(out, fus_faa)

cassette <- do.call(rbind, lapply(c("Af_AFUA_2G01220", "Af_AFUA_1G13300"), function(w) {
  d <- bl(fus_faa, file.path(SYN, paste0(w, ".faa")),
          file.path(SYN, paste0("rev_b4_fusion__vs__", w, ".tsv")))
  if (is.null(d)) return(NULL)
  wl <- windows[[w]]
  data.frame(fusion_protein = d$q, target_window = w,
             target_protein = d$s,
             target_locus = wl$locus_tag[match(d$s, wl$protein_id)],
             target_offset = wl$offset[match(d$s, wl$protein_id)],
             target_product = wl$product[match(d$s, wl$protein_id)],
             pident = d$pident, aln_len = d$len, evalue = d$evalue, bits = d$bits,
             stringsAsFactors = FALSE)
}))
write.csv(cassette, file.path(TAB, "rev_b4_rib1_classII_cassette_fusion.csv"), row.names = FALSE)
cat("\n== class-II two-gene cassette / fusion test ==\n")
print(cassette[order(cassette$fusion_protein, -cassette$bits), ], row.names = FALSE)
