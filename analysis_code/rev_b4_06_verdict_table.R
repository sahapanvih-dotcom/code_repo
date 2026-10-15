#!/usr/bin/env Rscript
# rev_b4_06_verdict_table.R
# Revision analysis B4, step 6: consolidate every independent line of evidence bearing on
# "the published Candida RIB1 validation corresponds to AFUA_1G13300, not AFUA_2G01220"
# into one auditable table, each row carrying the number it rests on and its source file.
BASE <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
OG   <- file.path(BASE, "analysis/orthogroup"); TAB <- file.path(BASE, "output/table")

# --- AU tests, parsed out of the IQ-TREE report files ---------------------------------
parse_au <- function(f, labels) {
  if (!file.exists(f)) return(NULL)
  l <- readLines(f); i <- grep("^Tree      logL", l)
  if (!length(i)) return(NULL)
  blk <- l[(i + 2):(i + 1 + length(labels))]
  num <- lapply(strsplit(trimws(gsub("[+-] ", " ", blk)), " +"), as.numeric)
  data.frame(tree = labels,
             logL   = vapply(num, `[`, numeric(1), 2),
             deltaL = vapply(num, `[`, numeric(1), 3),
             p_AU   = vapply(num, function(x) x[length(x)], numeric(1)),
             stringsAsFactors = FALSE)
}
labs <- c("unconstrained ML", "H1 constrained: Candida RIB1 + AFUA_1G13300 class",
          "H2 constrained: Candida RIB1 + AFUA_2G01220 class")
au_full <- parse_au(file.path(OG, "tree/rev_b4_autest.iqtree"), labs)
au_dom  <- parse_au(file.path(OG, "tree/rev_b4_dom_autest.iqtree"), labs)
au <- rbind(cbind(alignment = "full-length (mafft L-INS-i, trimal -automated1, 359 cols)", au_full),
            cbind(alignment = "PF00925 domain only (trimal -automated1, 213 cols)",      au_dom))
write.csv(au, file.path(TAB, "rev_b4_rib1_AU_topology_tests.csv"), row.names = FALSE)
cat("\n== AU tests ==\n"); print(au, row.names = FALSE)

top <- read.csv(file.path(TAB, "rev_b4_rib1_topology_tests.csv"), stringsAsFactors = FALSE)
syn <- read.csv(file.path(TAB, "rev_b4_rib1_synteny_within_lineage.csv"), stringsAsFactors = FALSE)
mem <- read.csv(file.path(TAB, "rev_b4_rib1_family_membership.csv"), stringsAsFactors = FALSE)
g1 <- mem[mem$acc == "XP_077660111.1", ]; g2 <- mem[mem$acc == "XP_749258.2", ]

ev <- rbind(
 data.frame(evidence = "Family completeness",
   result = sprintf("%d PF00925 proteins from %d proteomes; genome-wide hmmsearch of PF00925 over all 14 proteomes returned the identical set",
                    nrow(mem), length(unique(mem$species))),
   bears_on = "neither", source = "analysis/orthogroup/pfam/rev_b4_PF00925_vs_searchdb.tbl"),
 data.frame(evidence = "Copy number",
   result = "every Aspergillus genome carries both classes; C. albicans, C. tropicalis, S. cerevisiae and E. coli each carry exactly ONE PF00925 protein, always class I",
   bears_on = "AFUA_1G13300", source = "output/table/rev_b4_rib1_family_membership.csv"),
 data.frame(evidence = "Pfam domain architecture",
   result = sprintf("AFUA_1G13300 = %s (PF00925 seq E=%s), identical to Ct/Ca/Sc RIB1 and E. coli ribA; AFUA_2G01220 = %s (PF00925 seq E=%s, ~27 orders of magnitude weaker)",
                    g1$domain_architecture, format(g1$PF00925_seq_evalue, scientific = TRUE),
                    g2$domain_architecture, format(g2$PF00925_seq_evalue, scientific = TRUE)),
   bears_on = "AFUA_1G13300", source = "output/figure/rev_b4_rib1_gene_tree.pdf p5"),
 data.frame(evidence = "blastp identity vs CTRG_00799 (published Candida RIB1)",
   result = sprintf("AFUA_1G13300 %.1f%% id, bits %.0f, E=%s; AFUA_2G01220 %.1f%% id, bits %.0f, E=%s (%.1fx lower bitscore)",
                    g1$pident_vs_CtRIB1, g1$bitscore_vs_CtRIB1, format(g1$evalue_vs_CtRIB1, scientific = TRUE),
                    g2$pident_vs_CtRIB1, g2$bitscore_vs_CtRIB1, format(g2$evalue_vs_CtRIB1, scientific = TRUE),
                    g1$bitscore_vs_CtRIB1 / g2$bitscore_vs_CtRIB1),
   bears_on = "AFUA_1G13300", source = "analysis/orthogroup/blast/rev_b4_seeds_vs_searchdb.tsv"),
 data.frame(evidence = "Unrooted quartet {1G13300, 2G01220, CTRG_00799, ribA}",
   result = sprintf("%d/%d ML analyses = Q3 (co-orthologues, Aspergillus-lineage duplication); %d/%d = Q1 (1G13300 + Candida RIB1); 0/%d = Q2 (2G01220 + Candida RIB1). Identical under all 5 rootings.",
                    sum(grepl("Q3", top$quartet_vs_CtRIB1)), nrow(top),
                    sum(grepl("Q1", top$quartet_vs_CtRIB1)), nrow(top), nrow(top)),
   bears_on = "AFUA_1G13300 or co-orthology; never AFUA_2G01220",
   source = "output/table/rev_b4_rib1_topology_tests.csv"),
 data.frame(evidence = "AU topology test (full-length alignment)",
   result = sprintf("H1 p-AU=%.3f, H2 p-AU=%.3f, dLnL=%.3f - statistically indistinguishable, neither rejected",
                    au_full$p_AU[2], au_full$p_AU[3], abs(au_full$logL[2] - au_full$logL[3])),
   bears_on = "CANNOT RESOLVE", source = "output/table/rev_b4_rib1_AU_topology_tests.csv"),
 data.frame(evidence = "AU topology test (PF00925 domain alignment)",
   result = sprintf("H1 is the ML topology (dLnL=%.4f, p-AU=%.3f); H2 dLnL=%.2f, p-AU=%.2f - disfavoured but NOT rejected at 5%%",
                    au_dom$deltaL[2], au_dom$p_AU[2], au_dom$deltaL[3], au_dom$p_AU[3]),
   bears_on = "AFUA_1G13300 (weakly)", source = "output/table/rev_b4_rib1_AU_topology_tests.csv"),
 data.frame(evidence = "Branch support for (class I + Saccharomycotina RIB1) excluding class II",
   result = sprintf("recovered only in the PF00925-domain analysis, SH-aLRT/UFboot = %s (below the 80/95 thresholds)",
                    top$MRCA_1G13300_yeastRIB1_support[top$analysis == "domain_PF00925_only"]),
   bears_on = "AFUA_1G13300 (weak support)", source = "analysis/orthogroup/tree/rev_b4_domain.treefile"),
 data.frame(evidence = "Patristic distance to CTRG_00799",
   result = sprintf("AFUA_2G01220 is %.1fx-%.1fx further from the published Candida RIB1 than AFUA_1G13300 is, in every analysis",
                    min(top$patristic_ratio_2G_over_1G), max(top$patristic_ratio_2G_over_1G)),
   bears_on = "AFUA_1G13300", source = "output/table/rev_b4_rib1_topology_tests.csv"),
 data.frame(evidence = "Class-II stem branch length",
   result = "3.2-7.5 substitutions/site in the full-length trees: saturated, so class II is long-branch-attraction prone and its exact placement is NOT trustworthy",
   bears_on = "CAVEAT", source = "analysis/orthogroup/tree/*.treefile"),
 data.frame(evidence = "Synteny, Candida vs Candida around RIB1",
   result = sprintf("CTRG_00799 vs C. albicans RIB1: %d/10 flanking genes reciprocal (%.0f%%) - CTRG_00799 IS the locus the published C. albicans validation refers to",
                    syn$n_reciprocal_flanking_pairs[1], 100 * syn$flanking_synteny_frac[1]),
   bears_on = "validity of the Candida anchor", source = "output/table/rev_b4_rib1_synteny_within_lineage.csv"),
 data.frame(evidence = "Synteny, Aspergillus class I vs class I",
   result = sprintf("AFUA_1G13300 neighbourhood vs A. lentulus / A. nidulans / A. niger class I: %s flanking genes reciprocal of 10 - a deeply conserved Aspergillus locus",
                    paste(syn$n_reciprocal_flanking_pairs[syn$comparison_type == "Aspergillus class I vs class I"], collapse = "/")),
   bears_on = "AFUA_1G13300", source = "output/table/rev_b4_rib1_synteny_within_lineage.csv"),
 data.frame(evidence = "Synteny, class I vs class II in the same genome",
   result = "0/10 shared flanking genes - the two A. fumigatus loci are entirely distinct genomic positions, so neither is an annotation artefact of the other",
   bears_on = "neither", source = "output/table/rev_b4_rib1_synteny_within_lineage.csv"),
 data.frame(evidence = "Synteny, Aspergillus vs Candida/Saccharomyces",
   result = "0/10 shared flanking genes for BOTH A. fumigatus loci - microsynteny does not survive the Pezizomycotina/Saccharomycotina split, so synteny cannot adjudicate between the two paralogues directly",
   bears_on = "CANNOT RESOLVE", source = "output/table/rev_b4_rib1_synteny_summary.csv"),
 data.frame(evidence = "Class-II two-gene cassette",
   result = "A. terreus XP_001216775.1 and A. candidus XP_024676616.1 are single ORFs fusing AFUA_2G01220 (92.8%/87.4% id) to its immediate neighbour AFUA_2G01230, a DHDPS-family gene (80.5%/75.8% id); they hit AFUA_1G13300 at only ~32% over 88 aa. Class II is a Pezizomycotina-restricted gene module with no Saccharomycotina counterpart.",
   bears_on = "AFUA_1G13300", source = "output/table/rev_b4_rib1_classII_cassette_fusion.csv"),
 data.frame(evidence = "A. ochraceus tblastn (unannotated genome, independent check)",
   result = "C. tropicalis RIB1, C. albicans RIB1, S. cerevisiae RIB1 and E. coli ribA all hit the SAME A. ochraceus locus (SSOS01000005.1:2737740-2738963), which is also the best hit of AFUA_1G13300 (77.0% id); AFUA_2G01220 hits different contigs entirely",
   bears_on = "AFUA_1G13300", source = "analysis/orthogroup/blast/rev_b4_ochraceus_tblastn.tsv"),
 data.frame(evidence = "blastp -evalue reporting artefact (robustness of the family)",
   result = "blastp -evalue also sets the ungapped-extension cutoff, so a search run at -evalue 1e-5 can drop alignments whose final E is below 1e-5. Re-running the seed search at -evalue 10 and applying E<=1e-5 / qcov>=40 afterwards returns the IDENTICAL 29 proteins, and the genome-wide PF00925 hmmsearch (not subject to the artefact) returns the same 29. Under permissive reporting AFUA_1G13300 vs AFUA_2G01220 = 30.8% id over 91 aa (qcov 17-22%), E=9.4e-09, bits 45 - still ~8x weaker than AFUA_1G13300 vs CTRG_00799 (bits 350).",
   bears_on = "CAVEAT (checked; family unaffected)", source = "analysis/orthogroup/logs/rev_b4_versions_and_commands.txt"),
 data.frame(evidence = "Host (human) homology",
   result = "0 hits at E<=1e-3 for either A. fumigatus protein against the human proteome",
   bears_on = "neither", source = "analysis/orthogroup/blast/rev_b4_seeds_vs_human.tsv"),
 stringsAsFactors = FALSE)

ev$line <- seq_len(nrow(ev))
write.csv(ev[, c("line", "evidence", "result", "bears_on", "source")],
          file.path(TAB, "rev_b4_rib1_orthology_verdict.csv"), row.names = FALSE)
cat("\n== consolidated evidence ==\n")
for (i in seq_len(nrow(ev))) cat(sprintf("%2d. [%s] %s\n    %s\n", ev$line[i], ev$bears_on[i], ev$evidence[i], ev$result[i]))
