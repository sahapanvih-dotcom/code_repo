#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# Advisor round 3, task 4 (reproducibility): generate output/table/rev_c2_candidate_set_union.csv.
#
# The 26-gene lead-candidate set is the union of two lists. This script existed only as an
# interactive step until now, which left the file without a generator; the round-3 re-execution
# check flagged that, so the step is recorded here. Nothing in the file changes.
#
#   (1) The 18 Tier-1 candidates. DERIVED: every row of output/table/ms_table3_tier1.csv, with
#       that table's own pipeline rank and product string.
#
#   (2) The 12 genes carried in main Table 1. CURATED, and declared as such: this is the editorial
#       selection of genes the manuscript discusses by name, not the output of a filter. Four of
#       the 12 are also Tier-1 (AFUA_2G01220, AFUA_3G10310, AFUA_5G09750, AFUA_6G03540) and are
#       taken from (1); the remaining 8 are listed below in the order in which they entered the
#       table. That order is editorial. It is fixed here only so that `rank` is a stable sort key
#       for genes that have no pipeline rank, because they sit outside the ranked 815; rev_c2
#       re-sorts its output on triage and constraint percentile, so the value never reaches a
#       published number. Product strings are the Af293 annotations used throughout the project.
#
# Run: Rscript analysis_code/rev_c2_00_candidate_set_union.R  [out.csv]
#      (writes output/table/rev_c2_candidate_set_union.csv unless a path is given)
# ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table)})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
TAB <- "output/table"
args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else file.path(TAB, "rev_c2_candidate_set_union.csv")

tier1 <- fread(file.path(TAB, "ms_table3_tier1.csv"))[, .(locus_tag, product, rank, set = "tier1")]
setorder(tier1, rank)

# CURATED: the main Table 1 genes that are not Tier-1, in table-entry order (see header).
main_only <- data.table(
  locus_tag = c("AFUA_1G09840", "AFUA_1G13300", "AFUA_8G05210", "AFUA_2G03560",
                "AFUA_2G05060", "AFUA_1G02590", "AFUA_4G13510", "AFUA_4G10460"),
  product   = c("protein translocase subunit TIM44",
                "GTP cyclohydrolase II",
                "class I SAM-dependent methyltransferase",
                "putative two-component osmosensing histidine kinase (Bos1)",
                "alternative oxidase",
                "rDNA-binding RNA polymerase I transcriptional factor",
                "isocitrate lyase 1",
                "homocitrate synthase"),
  set       = "main_table1_only")
main_only[, rank := 9018L + .I]   # sort key only; these genes have no pipeline rank

stopifnot(!any(main_only$locus_tag %in% tier1$locus_tag))
out <- rbind(tier1, main_only)[, .(locus_tag, product, rank, set)]
stopifnot(nrow(out) == 26L, !anyDuplicated(out$locus_tag))
fwrite(out, OUT)
cat(sprintf("wrote %s: %d genes (%d tier1 + %d main_table1_only)\n",
            OUT, nrow(out), sum(out$set == "tier1"), sum(out$set == "main_table1_only")))
