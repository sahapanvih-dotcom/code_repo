#!/usr/bin/env Rscript
# rev_b4_02b_constraints.R
# Build the two competing topological hypotheses as IQ-TREE constraint trees.
#
#   H1 (manuscript's reassignment): the Saccharomycotina RIB1 proteins (C. tropicalis
#       CTRG_00799, C. albicans RIB1, S. cerevisiae RIB1) are sister to the Aspergillus
#       clade that contains AFUA_1G13300.
#   H2 (alternative):              they are sister to the Aspergillus clade that
#       contains AFUA_2G01220.
#
# Constraints are built from the two unambiguous domain-architecture classes only
# (bare PF00925 vs PF12471+PF00925). The two long DHDPS-fusion proteins and the short
# A. tubingensis copy are intentionally left OUT of the constraint so their placement
# stays free and cannot be used to manufacture either answer.
suppressPackageStartupMessages(library(ape))

BASE <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
OG   <- file.path(BASE, "analysis/orthogroup")

fam <- read.csv(file.path(OG, "rev_b4_family_coverage_domains.csv"), stringsAsFactors = FALSE)

OUTGROUP <- fam$tip[fam$tag == "Ecol"]
YEAST    <- fam$tip[fam$tag %in% c("Ctro", "Calb", "Scer")]
ASP      <- fam[fam$tag %in% c("Acal","Acan","Afla","Afum","Alen","Anid","Anig","Aory","Ater","Atub"), ]

# unambiguous classes: exactly 2 domains (GTP_CH_N+GTP_cyclohydro2) vs exactly PF00925,
# and length within the modal range of the class (excludes the 785/787 aa fusions and
# the 317 aa short copy)
clsII <- ASP$tip[ASP$domain_architecture == "GTP_CH_N+GTP_cyclohydro2" & ASP$aa_len < 600]
clsI  <- ASP$tip[ASP$domain_architecture == "GTP_cyclohydro2"          & ASP$aa_len > 350]

stopifnot(length(OUTGROUP) == 1, length(YEAST) == 3,
          any(grepl("XP_?077660111|XP077660111", clsI)),
          any(grepl("XP_?7492582|XP7492582",    clsII)))

grp <- function(x) paste0("(", paste(x, collapse = ","), ")")

# H1: ((classI, yeast), classII)  rooted on the bacterial outgroup
h1 <- sprintf("(%s,(%s,%s),%s);", OUTGROUP, grp(clsI), grp(YEAST), grp(clsII))
# H2: ((classII, yeast), classI)
h2 <- sprintf("(%s,(%s,%s),%s);", OUTGROUP, grp(clsII), grp(YEAST), grp(clsI))

writeLines(h1, file.path(OG, "tree/rev_b4_constraint_H1.tre"))
writeLines(h2, file.path(OG, "tree/rev_b4_constraint_H2.tre"))

info <- data.frame(
  hypothesis = c("H1", "H2"),
  statement = c("Saccharomycotina RIB1 sister to the AFUA_1G13300-type clade (manuscript)",
                "Saccharomycotina RIB1 sister to the AFUA_2G01220-type clade (alternative)"),
  n_classI_constrained  = length(clsI),
  n_classII_constrained = length(clsII),
  n_yeast = length(YEAST),
  n_free  = nrow(fam) - length(clsI) - length(clsII) - length(YEAST) - 1,
  stringsAsFactors = FALSE)
write.csv(info, file.path(OG, "tree/rev_b4_constraint_definitions.csv"), row.names = FALSE)

cat("classI (bare PF00925, Aspergillus) n =", length(clsI), "\n")
cat("classII (PF12471+PF00925, Aspergillus) n =", length(clsII), "\n")
cat("free (unconstrained) tips n =", info$n_free[1], "\n")
cat("H1:", h1, "\nH2:", h2, "\n")
