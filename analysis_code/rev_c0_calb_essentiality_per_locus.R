#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# Advisor round 2, WP5/WP6 support: one per-A.-fumigatus-locus summary of published
# C. albicans orthologue essentiality, so that every candidate-level statement draws on the
# same numbers as the benchmark in rev_b1_essentiality_benchmark.R.
#
# The identifier crosswalk and the per-source admissibility rules are reproduced verbatim
# from rev_b1_essentiality_benchmark.R sections 5 (they are documented there), because the
# advisor's gate requires one source of truth rather than a second, independently written
# mapping that could drift. Both collapse rules are emitted:
#   strict  within-source conflicts excluded (the benchmark's primary rule)
#   anyess  any ESSENTIAL row wins (the rule the earlier manuscript tables used)
# The output is cross-checked against the published ms_table10_intersection_essentiality.csv.
#
# IMPORTANT: this is published essentiality of the C. albicans ORTHOLOGUE. It is not, and
# must never be reported as, experimental validation of an A. fumigatus target.
#
# Run: Rscript analysis_code/rev_c0_calb_essentiality_per_locus.R
# ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table)})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
TAB <- "output/table"
say <- function(...) cat(sprintf(...), "\n", sep = "")

om  <- fread("analysis/crossgenus/ortho_map.tsv")
ess <- fread(cmd = "zcat analysis/crossgenus/calb_essentiality.tsv.gz")
idm <- fread("analysis/crossgenus/calb_idmap.tsv")

## ---- identifier crosswalk (verbatim from rev_b1 section 5) --------------------
setnames(idm, c("acc", "oln", "primary", "syn"))
idl  <- idm[, .(id = unlist(strsplit(oln, " "))), by = acc][id != ""]
norm <- function(x) gsub("_", "", sub("^CAALFM_", "", toupper(trimws(x))))
ess[, gid := norm(gene_id)]; idl[, gid := norm(id)]
cw <- ess[grepl("orf19=orf19", evidence) & grepl("^C[0-9]_", gene_id),
          .(gid_orf19 = norm(sub(".*orf19=(orf19[.0-9]+).*", "\\1", evidence)), gid = norm(gene_id))]
cw    <- unique(cw[gid_orf19 != "" & gid != ""])
extra <- unique(merge(cw, idl[, .(acc, gid)], by = "gid")[, .(acc, gid = gid_orf19)])
idl2  <- unique(rbind(idl[, .(acc, gid)], extra))
e <- merge(ess, idl2, by = "gid", allow.cartesian = TRUE)
say("essentiality rows mapped to a UniProt accession: %d of %d", nrow(e), nrow(ess))

## ---- per-source admissibility (verbatim from rev_b1) -------------------------
e[, het       := source == "CGD_phenotype" & grepl("heterozygous", evidence, ignore.case = TRUE)]
e[, redundant := source == "CGD_phenotype" & grepl("PMID: *(30377286|34764269)", evidence)]
src_defs <- list(
  CGD_phenotype      = quote(source == "CGD_phenotype" & !het),
  GRACEv1_OMeara2015 = quote(source == "GRACEv1_OMeara2015"),
  GRACEv1_Fu2021     = quote(source == "GRACEv1_Fu2021"),
  GRACEv2_Fu2021     = quote(source == "GRACEv2_Fu2021"),
  GRACEv3_Fu2025     = quote(source == "GRACEv3_Fu2025"),
  Segal2018          = quote(source == "Segal2018"))

per <- rbindlist(lapply(names(src_defs), function(sn) {
  sub <- e[eval(src_defs[[sn]])][call %in% c("ESSENTIAL", "NONESSENTIAL")]
  gl  <- sub[, .(nE = sum(call == "ESSENTIAL"), nN = sum(call == "NONESSENTIAL")), by = acc]
  gl[, `:=`(source = sn, conflict = nE > 0 & nN > 0,
            call_anyess = as.integer(nE > 0),
            call_strict = fifelse(nE > 0 & nN > 0, NA_integer_, as.integer(nE > 0)))]
  gl[, .(acc, source, conflict, call_anyess, call_strict)]
}))

summ <- per[, .(
  n_src_reporting_anyess = .N,
  n_src_essential_anyess = sum(call_anyess),
  n_src_reporting_strict = sum(!is.na(call_strict)),
  n_src_essential_strict = sum(call_strict, na.rm = TRUE),
  sources_essential      = paste(sort(source[call_anyess == 1]), collapse = ";")), by = acc]

## ---- attach to A. fumigatus loci --------------------------------------------
link <- om[!is.na(calb_acc) & calb_acc != "",
           .(af_locus, af_protein, calb_acc, calb_pident, calb_qcov)]
out <- merge(link, summ, by.x = "calb_acc", by.y = "acc", all.x = TRUE)
for (j in c("n_src_reporting_anyess", "n_src_essential_anyess",
            "n_src_reporting_strict", "n_src_essential_strict"))
  out[is.na(get(j)), (j) := 0L]
out[, evidence_class := "published C. albicans orthologue essentiality (NOT A. fumigatus validation)"]
setcolorder(out, c("af_locus", "af_protein", "calb_acc", "calb_pident", "calb_qcov"))
setorder(out, af_locus)
fwrite(out, file.path(TAB, "rev_c0_calb_essentiality_per_af_locus.csv"))
say("wrote rev_c0_calb_essentiality_per_af_locus.csv (%d A. fumigatus loci with a C. albicans orthologue)",
    nrow(out))

## ---- cross-check against the published table --------------------------------
pub <- fread(file.path(TAB, "ms_table10_intersection_essentiality.csv"))
chk <- merge(pub[, .(af_locus, pub_n_ess = n_ess, pub_n_call = n_call)],
             out[, .(af_locus, n_src_essential_anyess, n_src_reporting_anyess)],
             by = "af_locus")
chk[, agree := pub_n_ess == n_src_essential_anyess & pub_n_call == n_src_reporting_anyess]
fwrite(chk, file.path(TAB, "rev_c0_essentiality_crosscheck.csv"))
say("\ncross-check against ms_table10_intersection_essentiality.csv: %d of %d loci agree",
    sum(chk$agree), nrow(chk))
if (any(!chk$agree)) print(chk[agree == FALSE])
