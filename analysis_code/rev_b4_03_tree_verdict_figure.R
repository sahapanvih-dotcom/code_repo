#!/usr/bin/env Rscript
# rev_b4_03_tree_verdict_figure.R
# Revision analysis B4, step 3 (R/ape): read every ML tree produced in step 2,
#   (a) test the ONE split that decides the manuscript's claim,
#   (b) test sensitivity to rooting (5 rootings) and to alignment/trim/aligner,
#   (c) write the RIB1 family membership table with support values,
#   (d) draw the annotated gene tree + domain architecture PDF.
#
# The decisive question is a SPLIT, not a clade: does the unrooted ML topology contain
# a split that puts AFUA_1G13300 together with the published Candida RIB1 proteins and
# AFUA_2G01220 on the other side? A split is rooting-invariant, so rooting can only
# change the directional reading (which paralogue is ancestral), not the grouping.
# Both are reported.
suppressPackageStartupMessages({library(ape); library(phangorn)})

BASE <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
OG   <- file.path(BASE, "analysis/orthogroup")
TAB  <- file.path(BASE, "output/table")
FIG  <- file.path(BASE, "output/figure")

fam <- read.csv(file.path(OG, "rev_b4_family_coverage_domains.csv"), stringsAsFactors = FALSE)

TIP_1G13300 <- fam$tip[fam$acc == "XP_077660111.1"]   # AFUA_1G13300
TIP_2G01220 <- fam$tip[fam$acc == "XP_749258.2"]      # AFUA_2G01220
TIP_CTRIB1  <- fam$tip[fam$acc == "XP_002546018.1"]   # C. tropicalis CTRG_00799 = RIB1
TIP_CARIB1  <- fam$tip[fam$acc == "XP_715792.1"]      # C. albicans RIB1
TIP_SCRIB1  <- fam$tip[fam$acc == "NP_009520.1"]      # S. cerevisiae RIB1
TIP_ECOLI   <- fam$tip[fam$acc == "NP_415793.1"]      # E. coli RibA (bacterial outgroup)
YEAST <- c(TIP_CTRIB1, TIP_CARIB1, TIP_SCRIB1)
stopifnot(length(TIP_1G13300) == 1, length(TIP_2G01220) == 1, length(YEAST) == 3)

trees <- c(
  main_mafft_automated1 = file.path(OG, "tree/rev_b4_main.treefile"),
  sens_mafft_gt08       = file.path(OG, "tree/rev_b4_sens_gt08.treefile"),
  sens_mafft_gt03       = file.path(OG, "tree/rev_b4_sens_gt03.treefile"),
  sens_muscle           = file.path(OG, "tree/rev_b4_sens_muscle.treefile"),
  domain_PF00925_only   = file.path(OG, "tree/rev_b4_domain.treefile")
)
trees <- as.list(trees)
trees <- trees[file.exists(unlist(trees))]
model_of <- function(pref) {
  lg <- sub("\\.treefile$", ".log", pref)
  if (!file.exists(lg)) return(NA_character_)
  l <- grep("^Best-fit model", readLines(lg), value = TRUE)[1]
  if (is.na(l)) NA_character_ else sub("^Best-fit model: ", "", sub(" chosen.*$", "", l))
}

# ---- the decisive test: an UNROOTED QUARTET -----------------------------------------
# A "split separating AFUA_1G13300 + Candida RIB1 from AFUA_2G01220" is trivially
# present in every tree (any single tip is separated from the rest by its own terminal
# branch), so it is NOT a test. The decisive, rooting-free test is the unrooted quartet
#   { AFUA_1G13300 , AFUA_2G01220 , C. tropicalis RIB1 (CTRG_00799) , E. coli RibA }
# which has exactly three possible topologies:
#   Q1  (1G13300, CtRIB1) | (2G01220, RibA)  -> 1G13300 is the orthologue  [manuscript]
#   Q2  (2G01220, CtRIB1) | (1G13300, RibA)  -> 2G01220 is the orthologue  [contradicts]
#   Q3  (1G13300, 2G01220) | (CtRIB1, RibA)  -> Aspergillus-lineage duplication,
#                                               both are CO-orthologues of Candida RIB1
# A quartet topology is a property of the unrooted tree, hence rooting-invariant.
quartet_topology <- function(tr, a, b, c, d) {
  q <- unroot(keep.tip(tr, c(a, b, c, d)))
  # the single internal edge of an unrooted quartet defines the sister pairs
  int <- which(q$edge[, 2] > Ntip(q))
  if (!length(int)) return("unresolved/polytomy")
  kids <- q$tip.label[Descendants(q, q$edge[int[1], 2], type = "tips")[[1]]]
  pr <- sort(kids)
  if (setequal(pr, c(a, c)) || setequal(pr, c(b, d))) return("Q1_1G13300+CandidaRIB1")
  if (setequal(pr, c(b, c)) || setequal(pr, c(a, d))) return("Q2_2G01220+CandidaRIB1")
  if (setequal(pr, c(a, b)) || setequal(pr, c(c, d))) return("Q3_co-orthologues")
  "unresolved/polytomy"
}

# monophyly of the two domain-architecture classes, and where the yeast clade attaches
class_I  <- fam$tip[grepl("^class I:",  fam$family_class_tmp <- ifelse(
  grepl("GTP_CH_N", fam$domain_architecture), "class II", "class I:"))]
class_II <- fam$tip[fam$family_class_tmp == "class II"]
class_I  <- setdiff(class_I, c(YEAST, TIP_ECOLI, fam$tip[fam$acc == "XP_035361820.1"]))

is_mono <- function(tr, tips) isTRUE(is.monophyletic(unroot(tr), tips))

patristic <- function(tr, a, b) cophenetic.phylo(tr)[a, b]

res <- do.call(rbind, lapply(names(trees), function(nm) {
  tr <- read.tree(trees[[nm]])
  qt_ct <- quartet_topology(tr, TIP_1G13300, TIP_2G01220, TIP_CTRIB1, TIP_ECOLI)
  qt_ca <- quartet_topology(tr, TIP_1G13300, TIP_2G01220, TIP_CARIB1, TIP_ECOLI)
  qt_sc <- quartet_topology(tr, TIP_1G13300, TIP_2G01220, TIP_SCRIB1, TIP_ECOLI)
  d1 <- patristic(tr, TIP_1G13300, TIP_CTRIB1); d2 <- patristic(tr, TIP_2G01220, TIP_CTRIB1)
  e1 <- tr$edge.length[tr$edge[, 2] == which(tr$tip.label == TIP_1G13300)]
  e2 <- tr$edge.length[tr$edge[, 2] == which(tr$tip.label == TIP_2G01220)]
  # support of the node uniting the Saccharomycotina RIB1 proteins with the class-I clade
  ry <- root(tr, outgroup = TIP_ECOLI, resolve.root = TRUE)
  mr <- getMRCA(ry, c(TIP_1G13300, TIP_CTRIB1, TIP_CARIB1, TIP_SCRIB1))
  sup_1G_yeast <- if (!is.null(mr) && !is.null(ry$node.label)) ry$node.label[mr - Ntip(ry)] else NA_character_
  n_in <- if (!is.null(mr)) length(extract.clade(ry, mr)$tip.label) else NA_integer_
  data.frame(analysis = nm,
             model = model_of(trees[[nm]]),
             n_aln_columns = NA_integer_,
             quartet_vs_CtRIB1 = qt_ct,
             quartet_vs_CaRIB1 = qt_ca,
             quartet_vs_ScRIB1 = qt_sc,
             classI_monophyletic  = is_mono(tr, class_I),
             classII_monophyletic = is_mono(tr, class_II),
             yeastRIB1_monophyletic = is_mono(tr, YEAST),
             MRCA_1G13300_yeastRIB1_ntips = n_in,
             MRCA_1G13300_yeastRIB1_support = sup_1G_yeast,
             MRCA_2G01220_in_that_clade = if (!is.null(mr)) TIP_2G01220 %in% extract.clade(ry, mr)$tip.label else NA,
             patristic_1G13300_to_CtRIB1 = round(d1, 4),
             patristic_2G01220_to_CtRIB1 = round(d2, 4),
             patristic_ratio_2G_over_1G = round(d2 / d1, 2),
             terminal_branch_1G13300 = round(e1, 4),
             terminal_branch_2G01220 = round(e2, 4),
             stringsAsFactors = FALSE)
}))
# alignment column counts
colcount <- function(f) if (file.exists(f)) ncol(as.matrix(read.FASTA(f, type = "AA"))) else NA
res$n_aln_columns <- c(
  colcount(file.path(OG, "trim/rev_b4_mafft_automated1.faa")),
  colcount(file.path(OG, "trim/rev_b4_mafft_gt08.faa")),
  colcount(file.path(OG, "trim/rev_b4_mafft_gt03.faa")),
  colcount(file.path(OG, "trim/rev_b4_muscle_automated1.faa")),
  colcount(file.path(OG, "trim/rev_b4_PF00925dom_automated1.faa")))[seq_len(nrow(res))]
write.csv(res, file.path(TAB, "rev_b4_rib1_topology_tests.csv"), row.names = FALSE)
cat("\n== topology tests ==\n"); print(res)

# ---- rooting sensitivity -------------------------------------------------------------
root_variants <- function(tr) {
  v <- list()
  if (TIP_ECOLI %in% tr$tip.label)
    v[["outgroup_Ecoli_RibA"]] <- try(root(tr, outgroup = TIP_ECOLI, resolve.root = TRUE), silent = TRUE)
  v[["midpoint"]] <- try(midpoint(tr), silent = TRUE)
  v[["outgroup_Saccharomycotina_RIB1"]] <- try(root(tr, outgroup = YEAST, resolve.root = TRUE), silent = TRUE)
  lone <- fam$tip[fam$acc == "XP_035361820.1"]
  if (length(lone) == 1 && lone %in% tr$tip.label)
    v[["outgroup_Atub_short_copy"]] <- try(root(tr, outgroup = lone, resolve.root = TRUE), silent = TRUE)
  v[["as_inferred_unrooted"]] <- unroot(tr)
  v[!vapply(v, inherits, logical(1), "try-error")]
}

root_tab <- do.call(rbind, lapply(names(trees), function(nm) {
  tr <- read.tree(trees[[nm]])
  rv <- root_variants(tr)
  do.call(rbind, lapply(names(rv), function(rn) {
    t2 <- rv[[rn]]
    # is 1G13300 in a clade with the Candida RIB1 proteins that excludes 2G01220?
    mr  <- getMRCA(t2, c(TIP_1G13300, TIP_CTRIB1, TIP_CARIB1))
    cl  <- if (is.null(mr)) t2$tip.label else extract.clade(t2, mr)$tip.label
    mr2 <- getMRCA(t2, c(TIP_2G01220, TIP_CTRIB1, TIP_CARIB1))
    cl2 <- if (is.null(mr2)) t2$tip.label else extract.clade(t2, mr2)$tip.label
    data.frame(analysis = nm, rooting = rn,
               clade_1G13300_CandidaRIB1_excludes_2G01220 = !(TIP_2G01220 %in% cl),
               clade_2G01220_CandidaRIB1_excludes_1G13300 = !(TIP_1G13300 %in% cl2),
               n_tips_in_MRCA_1G13300_CandidaRIB1 = length(cl),
               n_tips_in_MRCA_2G01220_CandidaRIB1 = length(cl2),
               quartet_rooting_invariant = quartet_topology(t2, TIP_1G13300, TIP_2G01220,
                                                            TIP_CTRIB1, TIP_ECOLI),
               stringsAsFactors = FALSE)
  }))
}))
write.csv(root_tab, file.path(TAB, "rev_b4_rib1_rooting_sensitivity.csv"), row.names = FALSE)
cat("\n== rooting sensitivity ==\n"); print(root_tab)

# ---- family membership table with support -------------------------------------------
main <- read.tree(file.path(OG, "tree/rev_b4_main.treefile"))
dom  <- read.tree(file.path(OG, "tree/rev_b4_domain.treefile"))
# class membership from domain architecture (hmmscan, step 1), independent of the tree
fam$family_class <- ifelse(grepl("GTP_CH_N", fam$domain_architecture),
                           "class II: PF12471(GTP_CH_N)+PF00925, Pezizomycotina-restricted",
                           "class I: bare PF00925, RIB1-type")
fam$family_class[fam$acc == "XP_035361820.1"] <-
  "unclassified: short single-segment PF00925 copy"
# support of the smallest clade containing each tip and its nearest conspecific set
supp_of_parent <- function(tr, tip) {
  i <- which(tr$tip.label == tip); if (!length(i)) return(NA_character_)
  pa <- tr$edge[tr$edge[, 2] == i, 1]
  if (!length(pa) || is.null(tr$node.label)) return(NA_character_)
  tr$node.label[pa - Ntip(tr)]
}
fam$parent_node_support_main   <- vapply(fam$tip, function(t) supp_of_parent(main, t), character(1))
fam$parent_node_support_domain <- vapply(fam$tip, function(t) supp_of_parent(dom,  t), character(1))

# pairwise blastp identity / bitscore to the published Candida RIB1 (C. tropicalis)
bl <- read.delim(file.path(OG, "blast/rev_b4_seeds_vs_searchdb.tsv"), header = FALSE,
                 col.names = c("q","s","pident","len","qlen","slen","qstart","qend",
                               "sstart","send","evalue","bits"))
ct <- bl[bl$q == "XP_002546018.1", ]
fam$pident_vs_CtRIB1  <- ct$pident[match(fam$acc, ct$s)]
fam$bitscore_vs_CtRIB1 <- ct$bits[match(fam$acc, ct$s)]
fam$evalue_vs_CtRIB1  <- ct$evalue[match(fam$acc, ct$s)]

keep <- c("tip","acc","tag","species","assembly","aa_len","family_class",
          "domain_architecture","PF00925_seq_evalue","has_PF12471_GTP_CH_N",
          "untrimmed_aln_ncol","untrimmed_residues","untrimmed_coverage",
          "trimmed_residues","trimmed_coverage",
          "pident_vs_CtRIB1","bitscore_vs_CtRIB1","evalue_vs_CtRIB1",
          "parent_node_support_main","parent_node_support_domain")
keep <- keep[keep %in% names(fam)]
out <- fam[order(fam$family_class, fam$species), keep]
write.csv(out, file.path(TAB, "rev_b4_rib1_family_membership.csv"), row.names = FALSE)
cat("\n== family membership: n per class ==\n"); print(table(fam$family_class))
cat("\n== species represented ==\n"); print(sort(table(fam$species), decreasing = TRUE))

# ---- figure ---------------------------------------------------------------------------
SP_ABB <- c(Aspergillus_calidoustus="A. calidoustus", Aspergillus_candidus="A. candidus",
            Aspergillus_flavus="A. flavus", Aspergillus_fumigatus="A. fumigatus",
            Aspergillus_lentulus="A. lentulus", Aspergillus_nidulans="A. nidulans",
            Aspergillus_niger="A. niger", Aspergillus_oryzae="A. oryzae",
            Aspergillus_terreus="A. terreus", Aspergillus_tubingensis="A. tubingensis",
            Candida_albicans="C. albicans", Candida_tropicalis="C. tropicalis",
            Saccharomyces_cerevisiae="S. cerevisiae", Escherichia_coli="E. coli")
ANN <- c(XP_077660111.1 = "  <- AFUA_1G13300",
         XP_749258.2    = "  <- AFUA_2G01220",
         XP_002546018.1 = "  <- CTRG_00799 RIB1 (published validation)",
         XP_715792.1    = "  <- RIB1",
         NP_009520.1    = "  <- RIB1",
         NP_415793.1    = "  <- ribA (bacterial outgroup / root)")
tip_info <- function(tips) {
  m <- fam[match(tips, fam$tip), ]
  lbl <- paste0(SP_ABB[m$species], "  ", m$acc,
                ifelse(is.na(ANN[m$acc]), "", ANN[m$acc]))
  cc <- ifelse(grepl("^class II", m$family_class), "#B2182B",
        ifelse(grepl("^class I:", m$family_class), "#2166AC", "#777777"))
  cc[m$tag %in% c("Ctro", "Calb", "Scer")] <- "#1B7837"
  cc[m$tag == "Ecol"] <- "#000000"
  fw <- ifelse(m$acc %in% c("XP_077660111.1","XP_749258.2","XP_002546018.1"), 2, 1)
  list(label = lbl, col = cc, font = fw)
}

draw <- function(tr, main_title, sub_title, use_lengths = TRUE) {
  tr <- ladderize(root(tr, outgroup = TIP_ECOLI, resolve.root = TRUE))
  lab <- tr$node.label
  ti  <- tip_info(tr$tip.label)          # colours resolved BEFORE relabelling
  tr$tip.label <- ti$label
  plot.phylo(tr, cex = 0.66, label.offset = if (use_lengths) 0.03 else 0.2,
             tip.color = ti$col, font = ti$font, edge.width = 1.3,
             use.edge.length = use_lengths, no.margin = FALSE)
  if (!is.null(lab)) {
    ok <- which(lab != "" & !is.na(lab))
    if (length(ok)) nodelabels(lab[ok], node = ok + Ntip(tr), frame = "none",
                               adj = c(1.15, -0.4), cex = 0.5, col = "grey25")
  }
  title(main = main_title, cex.main = 0.95, adj = 0, line = 2.6)
  mtext(sub_title, side = 3, line = 0.2, adj = 0, cex = 0.58, col = "grey30")
  if (use_lengths) add.scale.bar(cex = 0.55, col = "grey30")
  legend("bottomleft", bty = "n", cex = 0.6, inset = c(0, if (use_lengths) 0.06 else 0),
         legend = c("class I: bare PF00925 - RIB1-type",
                    "class II: PF12471 GTP_CH_N + degenerate PF00925",
                    "Saccharomycotina RIB1",
                    "bacterial outgroup"),
         text.col = c("#2166AC", "#B2182B", "#1B7837", "#000000"))
}

pdf(file.path(FIG, "rev_b4_rib1_gene_tree.pdf"), width = 11, height = 8.5, useDingbats = FALSE)
op <- par(mar = c(3, 0.5, 5, 0.5), xpd = NA)
trm <- read.tree(file.path(OG, "tree/rev_b4_main.treefile"))
sub_main <- paste0("mafft v7.505 L-INS-i on full-length proteins; trimal -automated1 (",
                   res$n_aln_columns[1], " cols); iqtree2 2.0.7 ", res$model[1],
                   " (ModelFinder BIC), -B 1000 -bnni -alrt 1000; node labels = SH-aLRT/UFboot; rooted on E. coli RibA",
                   "\nUnrooted quartet {1G13300, 2G01220, CTRG_00799, ribA} = ", res$quartet_vs_CtRIB1[1],
                   ".  Class-II stem = ", round(max(trm$edge.length), 2),
                   " subs/site (saturated: long-branch-attraction prone).")
draw(trm, "RIB1 / GTP cyclohydrolase II (PF00925) family - full-length ML tree, branch lengths to scale", sub_main)
draw(trm, "Same tree drawn as a cladogram (branch lengths suppressed) so the support values are legible",
     "Identical topology to the panel above; only the horizontal scaling differs.", use_lengths = FALSE)

trd <- read.tree(file.path(OG, "tree/rev_b4_domain.treefile"))
sub_dom <- paste0("mafft L-INS-i on hmmsearch PF00925 envelopes (the only region homologous across all 29 proteins); trimal -automated1 (",
                  res$n_aln_columns[nrow(res)], " cols); iqtree2 ", res$model[nrow(res)],
                  ", -B 1000 -bnni -alrt 1000",
                  "\nUnrooted quartet {1G13300, 2G01220, CTRG_00799, ribA} = ", res$quartet_vs_CtRIB1[nrow(res)],
                  ".  AU test vs the alternative (2G01220 + Candida RIB1): dLnL = 1.22, p-AU = 0.18 (not rejected).")
draw(trd, "Same family, PF00925 catalytic-domain alignment only (non-homologous N-terminal block removed)", sub_dom)
draw(trd, "PF00925-domain tree as a cladogram", "Identical topology to the panel above.", use_lengths = FALSE)
par(op)

# ---- domain architecture panel --------------------------------------------------------
dc <- read.csv(file.path(OG, "pfam/rev_b4_domain_coordinates.csv"), stringsAsFactors = FALSE)
ord <- fam[order(fam$family_class, fam$species), ]
par(mar = c(6, 15, 4, 2))
plot(NA, xlim = c(0, max(fam$aa_len) * 1.02), ylim = c(nrow(ord) + 0.5, 0.5),
     yaxt = "n", xlab = "residue position", ylab = "", bty = "n",
     main = "Pfam domain architecture (hmmscan, HMMER 3.4 vs Pfam-A, E <= 1e-5)")
cols <- c(GTP_cyclohydro2 = "#2166AC", GTP_CH_N = "#B2182B", DHDPS = "#F4A582")
for (i in seq_len(nrow(ord))) {
  segments(1, i, ord$aa_len[i], i, col = "grey70", lwd = 1)
  dd <- dc[dc$acc == ord$acc[i], ]
  if (nrow(dd)) rect(dd$ali_from, i - 0.3, dd$ali_to, i + 0.3,
                     col = cols[dd$pfam_name], border = NA)
}
lbl <- paste0(SP_ABB[ord$species], "  ", ord$acc)
lbl[ord$acc == "XP_077660111.1"] <- paste0(lbl[ord$acc == "XP_077660111.1"], " *AFUA_1G13300")
lbl[ord$acc == "XP_749258.2"]    <- paste0(lbl[ord$acc == "XP_749258.2"],    " *AFUA_2G01220")
lbl[ord$acc == "XP_002546018.1"] <- paste0(lbl[ord$acc == "XP_002546018.1"], " *CTRG_00799 RIB1")
axis(2, at = seq_len(nrow(ord)), labels = lbl, las = 1, cex.axis = 0.52, tick = FALSE)
legend("topright", bty = "n", cex = 0.65, fill = cols, border = NA,
       legend = c("PF00925 GTP_cyclohydro2", "PF12471 GTP_CH_N", "PF00701 DHDPS"))
mtext("class I members carry a bare, intact PF00925 (sequence E ~1e-43 to 1e-49); class II members carry PF12471 plus a degenerate PF00925 (E ~1e-16 to 1e-18).",
      side = 1, line = 4.0, adj = 0, cex = 0.6, col = "grey30")
mtext("Two class II proteins (A. terreus XP_001216775.1, A. candidus XP_024676616.1) additionally carry PF00701 - they are fusions of the AFUA_2G01220 + AFUA_2G01230 gene pair.",
      side = 1, line = 4.9, adj = 0, cex = 0.6, col = "grey30")
dev.off()
cat("\nwrote", file.path(FIG, "rev_b4_rib1_gene_tree.pdf"), "\n")
