#!/usr/bin/env Rscript
# =============================================================================
# rev_b4_05_tier1_orthogroup_check.R
# Revision analysis B4.05 - orthogroup-aware orthology confidence for the
# 18 Tier-1 candidates of output/table/ms_table3_tier1.csv
#
# Primary language: R (project convention). The sequence searches themselves are
# delegated to analysis_code/rev_b4_05a_tier1_searches.sh because the correct
# tools for them are blastp / hmmscan / seqkit; all parsing, flagging,
# statistics and table writing are done here in R.
#
# OUTPUTS
#   output/table/rev_b4_tier1_orthology_confidence.csv
#   output/table/rev_b4_tier1_paralog_detail.csv
# INTERMEDIATES / LOGS
#   analysis/orthogroup/tier1/{fasta,db,blast,pfam,logs}/
#   analysis/orthogroup/tier1/rev_b4_tier1_versions.txt
#
# -----------------------------------------------------------------------------
# EXPLICIT CONFIDENCE RULE (also written into the confidence_rule CSV column)
# -----------------------------------------------------------------------------
# Measurements per candidate (all from blastp E<=1e-5 / hmmscan E<=1e-5):
#   qcov            = BLAST qcovs/100 of the forward best hit (query coverage)
#   margin          = bitscore(best hit) / bitscore(second-best distinct subject);
#                     Inf when only one hit passes E<=1e-5 (unambiguous)
#   reciprocal      = blastp(best Candida hit -> af293) top-1 returns the same
#                     A. fumigatus accession
#   architecture    = ordered Pfam-A domain string (i-Evalue<=1e-5, overlapping
#                     domains resolved best-score-first), tandem repeats of the
#                     same domain collapsed before comparison
#
# RISK FLAGS (each counts once toward n_risk_flags):
#   F1  paralogy                 : accession in has_paralog.txt OR n_af_inparalogs>=1
#                                  (the two are merged so paralogy is not counted twice)
#   F2  ct_low_qcov              : C. tropicalis forward qcov < 0.5
#   F3  ct_ambiguous_best        : C. tropicalis margin < 1.5
#   F4  ct_non_reciprocal        : C. tropicalis reciprocity FALSE
#   F5  calb_low_qcov            : C. albicans forward qcov < 0.5
#   F6  calb_ambiguous_best      : C. albicans margin < 1.5
#   F7  calb_non_reciprocal      : C. albicans reciprocity FALSE
#   F8a domain_mismatch_disjoint : architectures differ AND share NO Pfam domain name
#   F8b domain_partial_mismatch  : architectures differ but share >=1 Pfam domain name
#                                  (domain gain/loss, or Pfam using different models of
#                                  the same repeat family). Pfam clan assignments are NOT
#                                  available locally (Pfam-A.hmm carries no CL lines and
#                                  there is no Pfam-A.clans.tsv), so clan-level
#                                  equivalence cannot be computed; the shared-domain-name
#                                  overlap (domain_jaccard) is used as the objective
#                                  substitute. Only F8a forces UNRESOLVED.
#   F9  pfam_uninformative       : neither protein has any Pfam-A domain at E<=1e-5
#                                  (concordance cannot be tested)
#   F10 no_ct_hit                : no C. tropicalis hit at E<=1e-5
#   F11 no_calb_hit              : no C. albicans hit at E<=1e-5
#   F12 orthomap_no_ortholog     : ortho_map.tsv row carries no C. tropicalis ortholog
#
# CALL (evaluated in this order):
#   UNRESOLVED : no hit in both Candida species, OR non-reciprocal in both,
#                OR disjoint domain architecture (F8a)
#   HIGH       : reciprocal in BOTH species, no paralogy, qcov>=0.5 in both,
#                margin>=1.5 in both, architectures match -> n_risk_flags == 0
#   MEDIUM     : reciprocal in both species AND exactly one risk flag
#   LOW        : everything else (>=2 risk flags, or non-reciprocal in one species)
#
# assignment_changed = recomputed C. tropicalis best-hit locus differs from the
#   ct_locus recorded in ortho_map.tsv (NA/NA counts as unchanged); the direction
#   is given in assignment_change_type, and assignment_changed_confident is TRUE only
#   when the newly assigned C. tropicalis protein is also the reciprocal best hit.
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

# -----------------------------------------------------------------------------
# NOTE ON THE SEARCH THRESHOLD (important, and the reason the searches are run
# permissively): blastp uses -evalue both as the reporting threshold AND to set
# the preliminary ungapped-extension score cutoff, so running blastp at
# -evalue 1e-5 silently discards alignments whose FINAL E-value is far below
# 1e-5. Verified on this data set: A. fumigatus XP_749067.1 (AFUA_7G04420) vs
# C. tropicalis XP_002548511.1 aligns with E = 1.45e-09 over 70% of the query,
# yet is reported only when blastp is run with -evalue >= 1e-3. The companion
# script therefore reports at -evalue 10 and the E <= 1e-5 criterion is applied
# here in R. hmmscan is handled the same way (reported at the HMMER default,
# filtered at domain i-Evalue <= 1e-5 here).
# -----------------------------------------------------------------------------
EVAL_MAX <- 1e-5

BASE   <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
# The candidate set, the working directory and the output suffix are overridable so the
# same audit can be re-run over a wider set (advisor round 2, WP2) without touching the
# published Tier-1 outputs. Defaults reproduce the original run exactly.
CAND_CSV <- Sys.getenv("REV_B4_CAND_CSV", unset = "")
SUFFIX   <- Sys.getenv("REV_B4_OUT_SUFFIX", unset = "")
WD     <- file.path(BASE, Sys.getenv("REV_B4_WORKDIR", unset = "analysis/orthogroup/tier1"))
FA     <- file.path(WD, "fasta"); BL <- file.path(WD, "blast")
PF     <- file.path(WD, "pfam");  LG <- file.path(WD, "logs")
OUTTAB <- file.path(BASE, "output/table")
for (d in c(FA, BL, PF, LG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

F_TIER1   <- if (nzchar(CAND_CSV)) CAND_CSV else file.path(BASE, "output/table/ms_table3_tier1.csv")
F_ORTHO   <- file.path(BASE, "analysis/crossgenus/ortho_map.tsv")
F_AFFAA   <- file.path(BASE, "analysis/target_discovery/hom/af293.faa")
F_PARA    <- file.path(BASE, "analysis/target_discovery/hom/has_paralog.txt")
F_CTFAA   <- file.path(BASE, "analysis/candida/ref/protein.faa")
F_CTGFF   <- file.path(BASE, "analysis/candida/ref/genomic.gff")
F_CALBFAA <- file.path(BASE, "analysis/orthogroup/ref/calb_sc5314.faa")
F_CALBMAP <- file.path(BASE, "analysis/crossgenus/calb_idmap.tsv")
F_SEARCH  <- file.path(BASE, "analysis_code/rev_b4_05a_tier1_searches.sh")

stopifnot(all(file.exists(c(F_TIER1, F_ORTHO, F_AFFAA, F_PARA, F_CTFAA,
                            F_CTGFF, F_CALBFAA, F_CALBMAP, F_SEARCH))))

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

CONF_RULE <- paste(
  "UNRESOLVED if (no Candida hit in both species) or (non-reciprocal in both) or disjoint Pfam architecture (no domain name in common);",
  "HIGH if reciprocal in both C. tropicalis and C. albicans, no paralogy, forward qcov>=0.5 in both,",
  "bitscore margin to second-best >=1.5 in both, and Pfam architectures match (n_risk_flags==0);",
  "MEDIUM if reciprocal in both species and exactly 1 risk flag;",
  "LOW otherwise (>=2 risk flags, or non-reciprocal in one species).",
  "A differing architecture that still shares >=1 Pfam domain counts as one risk flag (domain_partial_mismatch) rather than forcing UNRESOLVED,",
  "because Pfam clan assignments are not available locally and repeat families are modelled inconsistently between the two proteins.",
  "All blastp/hmmscan significance thresholds E<=1e-5, applied in R (searches are reported permissively",
  "because blastp -evalue also sets the extension cutoff and would lose true hits);",
  "in-paralog = distinct af293 protein with E<=1e-5 and query coverage (qcovs) >=50.")

# ---------------------------------------------------------------- 1. inputs ----
msg("reading Tier-1 table and ortho_map")
tier1 <- read.csv(F_TIER1)
N_EXPECT <- if (nzchar(CAND_CSV)) nrow(tier1) else 18L
stopifnot(nrow(tier1) == N_EXPECT, all(c("rank", "locus_tag", "product") %in% names(tier1)))

om <- read.delim(F_ORTHO, header = TRUE, quote = "", colClasses = "character",
                 na.strings = c("", "NA", "-"))
stopifnot(all(c("af_protein", "af_locus", "ct_protein", "ct_locus",
                "ct_pident", "ct_qcov", "calb_acc", "calb_pident", "calb_qcov") %in% names(om)))

cand <- merge(tier1[, c("rank", "locus_tag", "product")], om,
              by.x = "locus_tag", by.y = "af_locus", all.x = TRUE, sort = FALSE)
cand <- cand[order(cand$rank), ]
stopifnot(nrow(cand) == N_EXPECT, !any(is.na(cand$af_protein)))
msg("resolved ", nrow(cand), " locus_tags to af_protein accessions")

writeLines(cand$af_protein, file.path(FA, "tier1_af_accessions.txt"))
system2("sync")  # project volume is NFS; make the list visible to the search script

# ----------------------------------------------- 2. run the external searches --
need <- c(file.path(BL, c("cand_vs_af293.tsv", "cand_vs_ctrop.tsv", "cand_vs_calb.tsv",
                          "ctrop_besthits_vs_af293.tsv", "calb_besthits_vs_af293.tsv")),
          file.path(PF, c("tier1_candidates_pfam.domtbl", "ctrop_besthits_pfam.domtbl")))
if (!all(file.exists(need)) || Sys.getenv("REV_B4_FORCE") == "1") {
  msg("launching rev_b4_05a_tier1_searches.sh (blastp / hmmscan / seqkit)")
  rc <- system2("bash", c(shQuote(F_SEARCH), shQuote(BASE)),
                stdout = file.path(LG, "rev_b4_05a_stdout.log"),
                stderr = file.path(LG, "rev_b4_05a_stderr.log"))
  if (rc != 0) stop("rev_b4_05a_tier1_searches.sh failed with status ", rc)
} else {
  msg("search outputs already present - reusing (set REV_B4_FORCE=1 to redo)")
}
# NFS metadata on this host can lag behind the writer process; wait for the files
for (k in 1:20) { if (all(file.exists(need))) break; Sys.sleep(1) }
stopifnot(all(file.exists(need)))

# ------------------------------------------------------------- 3. helpers -----
BCOLS <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore",
           "qlen", "slen", "qcovs")

read_blast <- function(f) {
  if (file.info(f)$size == 0) {
    x <- as.data.frame(setNames(replicate(length(BCOLS), character(0), simplify = FALSE), BCOLS))
    for (k in c("pident","length","evalue","bitscore","qlen","slen","qcovs")) x[[k]] <- numeric(0)
    return(x)
  }
  x <- read.delim(f, header = FALSE, col.names = BCOLS, quote = "", comment.char = "")
  # significance cut applied here, NOT by blastp (see the note at the top)
  x <- x[x$evalue <= EVAL_MAX, , drop = FALSE]
  # one row per (query, subject): the best-scoring HSP (qcovs is per subject pair)
  x <- x[order(x$qseqid, -x$bitscore), ]
  x[!duplicated(x[, c("qseqid", "sseqid")]), ]
}

fasta_products <- function(f) {
  h <- grep("^>", readLines(f, warn = FALSE), value = TRUE)
  acc <- sub("^>(\\S+).*$", "\\1", h)
  prod <- sub("^>\\S+\\s*", "", h)
  prod <- sub("\\s*\\[[^]]*\\]\\s*$", "", prod)
  setNames(prod, acc)
}

# ordered Pfam architecture from an hmmscan --domtblout file
parse_domtbl <- function(f, max_ievalue = EVAL_MAX) {
  ln <- readLines(f, warn = FALSE)
  ln <- ln[!grepl("^#", ln) & nzchar(trimws(ln))]
  if (!length(ln)) return(data.frame(query = character(0), domain = character(0),
                                     ali_from = numeric(0), score = numeric(0)))
  p <- strsplit(trimws(ln), "\\s+")
  d <- data.frame(domain   = vapply(p, `[`, "", 1),
                  query    = vapply(p, `[`, "", 4),
                  ievalue  = as.numeric(vapply(p, `[`, "", 13)),
                  score    = as.numeric(vapply(p, `[`, "", 14)),
                  ali_from = as.numeric(vapply(p, `[`, "", 18)),
                  ali_to   = as.numeric(vapply(p, `[`, "", 19)))
  d[d$ievalue <= max_ievalue, ]
}

# greedy best-score-first overlap resolution, then order along the query
architecture <- function(d, q) {
  dq <- d[d$query == q, ]
  if (!nrow(dq)) return("none")
  dq <- dq[order(-dq$score), ]
  keep <- dq[0, ]
  for (i in seq_len(nrow(dq))) {
    r <- dq[i, ]
    ok <- TRUE
    if (nrow(keep)) {
      ov <- pmin(r$ali_to, keep$ali_to) - pmax(r$ali_from, keep$ali_from) + 1
      ov[ov < 0] <- 0
      if (any(ov / (r$ali_to - r$ali_from + 1) > 0.5)) ok <- FALSE
    }
    if (ok) keep <- rbind(keep, r)
  }
  keep <- keep[order(keep$ali_from), ]
  paste(keep$domain, collapse = "-")
}
collapse_tandem <- function(s) {
  if (is.na(s) || s == "none") return(s)
  v <- strsplit(s, "-", fixed = TRUE)[[1]]
  paste(v[c(TRUE, v[-1] != v[-length(v)])], collapse = "-")
}

# ---------------------------------------------------------- 4. read searches ---
msg("parsing BLAST and hmmscan output")
b_af    <- read_blast(file.path(BL, "cand_vs_af293.tsv"))
b_ct    <- read_blast(file.path(BL, "cand_vs_ctrop.tsv"))
b_calb  <- read_blast(file.path(BL, "cand_vs_calb.tsv"))
r_ct    <- read_blast(file.path(BL, "ctrop_besthits_vs_af293.tsv"))
r_calb  <- read_blast(file.path(BL, "calb_besthits_vs_af293.tsv"))

af_prod   <- fasta_products(F_AFFAA)
ct_prod   <- fasta_products(F_CTFAA)
calb_prod <- fasta_products(F_CALBFAA)
paralog_list <- readLines(F_PARA, warn = FALSE)
paralog_list <- paralog_list[nzchar(paralog_list)]

# C. tropicalis protein_id -> locus_tag from the RefSeq GFF
msg("building C. tropicalis protein_id -> locus_tag map from genomic.gff")
gff <- readLines(F_CTGFF, warn = FALSE)
gff <- grep("\tCDS\t", gff, fixed = FALSE, value = TRUE)
att <- sub("^(?:[^\t]*\t){8}", "", gff)
pid <- sub("^.*protein_id=([^;]+).*$", "\\1", att); pid[!grepl("protein_id=", att)] <- NA
ltg <- sub("^.*locus_tag=([^;]+).*$", "\\1", att); ltg[!grepl("locus_tag=", att)] <- NA
ok  <- !is.na(pid) & !is.na(ltg)
ct_locus_of <- tapply(ltg[ok], pid[ok], function(z) z[1])
msg("  mapped ", length(ct_locus_of), " C. tropicalis protein accessions")

# C. albicans: RefSeq (local proteome) vs UniProt (ortho_map) are different ID
# spaces. Bridge as far as the local files allow: CAALFM_ locus from the RefSeq
# header where present, CAALFM_ locus from calb_idmap.tsv for the UniProt accession.
calb_locus_of_xp <- local({
  h <- grep("^>", readLines(F_CALBFAA, warn = FALSE), value = TRUE)
  acc <- sub("^>(\\S+).*$", "\\1", h)
  loc <- ifelse(grepl("CAALFM_[A-Za-z0-9]+", h),
                sub("^.*?(CAALFM_[A-Za-z0-9]+).*$", "\\1", h), NA_character_)
  setNames(loc, acc)
})
idmap <- read.delim(F_CALBMAP, header = TRUE, quote = "", colClasses = "character",
                    na.strings = c("", "NA"))
names(idmap)[1:3] <- c("entry", "ordered_locus", "gene_primary")
calb_locus_of_up <- setNames(sub("\\s.*$", "", idmap$ordered_locus), idmap$entry)

# ------------------------------------------------------ 5. per-candidate stats --
dom_af <- parse_domtbl(file.path(PF, "tier1_candidates_pfam.domtbl"))
dom_ct <- parse_domtbl(file.path(PF, "ctrop_besthits_pfam.domtbl"))

best_two <- function(tab, q) {
  s <- tab[tab$qseqid == q, ]
  if (!nrow(s)) return(NULL)
  s <- s[order(-s$bitscore), ]
  list(best = s[1, ], second = if (nrow(s) >= 2) s[2, ] else NULL, n = nrow(s))
}
recip_top <- function(tab, q) {
  s <- tab[tab$qseqid == q, ]
  if (!nrow(s)) return(NA_character_)
  s[which.max(s$bitscore), "sseqid"]
}

res <- list(); para_rows <- list()

for (i in seq_len(nrow(cand))) {
  cr <- cand[i, ]
  q  <- cr$af_protein

  ## --- paralogy ------------------------------------------------------------
  ip <- b_af[b_af$qseqid == q & b_af$sseqid != q & b_af$evalue <= EVAL_MAX & b_af$qcovs >= 50, ]
  ip <- ip[order(-ip$bitscore), ]
  n_ip <- nrow(ip)
  if (n_ip) {
    para_rows[[length(para_rows) + 1]] <- data.frame(
      locus_tag        = cr$locus_tag,
      af_protein       = q,
      inparalog_acc    = ip$sseqid,
      pident           = ip$pident,
      qcov             = ip$qcovs / 100,
      bitscore         = ip$bitscore,
      evalue           = ip$evalue,
      inparalog_product = unname(af_prod[ip$sseqid]))
  }
  has_par <- q %in% paralog_list

  ## --- forward / reciprocal: C. tropicalis ---------------------------------
  ct <- best_two(b_ct, q)
  if (is.null(ct)) {
    ct_acc <- NA_character_; ct_loc_new <- NA_character_; ct_pid <- NA_real_
    ct_qc <- NA_real_; ct_bs <- NA_real_; ct_mg <- NA_real_; ct_rec <- NA
    ct_n <- 0L; ct_back <- NA_character_
  } else {
    ct_acc <- ct$best$sseqid
    ct_loc_new <- unname(ct_locus_of[ct_acc]); if (is.null(ct_loc_new)) ct_loc_new <- NA_character_
    ct_pid <- ct$best$pident; ct_qc <- ct$best$qcovs / 100; ct_bs <- ct$best$bitscore
    ct_mg  <- if (is.null(ct$second)) Inf else ct$best$bitscore / ct$second$bitscore
    ct_n   <- ct$n
    ct_back <- recip_top(r_ct, ct_acc)
    ct_rec <- !is.na(ct_back) && ct_back == q
  }

  ## --- forward / reciprocal: C. albicans -----------------------------------
  cb <- best_two(b_calb, q)
  if (is.null(cb)) {
    cb_acc <- NA_character_; cb_pid <- NA_real_; cb_qc <- NA_real_; cb_bs <- NA_real_
    cb_mg <- NA_real_; cb_rec <- NA; cb_n <- 0L; cb_back <- NA_character_
  } else {
    cb_acc <- cb$best$sseqid
    cb_pid <- cb$best$pident; cb_qc <- cb$best$qcovs / 100; cb_bs <- cb$best$bitscore
    cb_mg  <- if (is.null(cb$second)) Inf else cb$best$bitscore / cb$second$bitscore
    cb_n   <- cb$n
    cb_back <- recip_top(r_calb, cb_acc)
    cb_rec <- !is.na(cb_back) && cb_back == q
  }

  ## --- Pfam architecture concordance ---------------------------------------
  arch_af <- architecture(dom_af, q)
  arch_ct <- if (is.na(ct_acc)) NA_character_ else architecture(dom_ct, ct_acc)
  a1 <- collapse_tandem(arch_af); a2 <- collapse_tandem(arch_ct)
  if (is.na(a2)) {
    dom_match <- NA; pfam_uninf <- FALSE
  } else if (a1 == "none" && a2 == "none") {
    dom_match <- NA; pfam_uninf <- TRUE
  } else {
    dom_match <- identical(a1, a2); pfam_uninf <- FALSE
  }
  dom_match_strict <- if (is.na(arch_ct)) NA else identical(arch_af, arch_ct)
  # domain-name set overlap: objective qualifier for a mismatch (no local clan data)
  s1 <- setdiff(unique(strsplit(a1, "-", fixed = TRUE)[[1]]), "none")
  s2 <- if (is.na(a2)) character(0) else setdiff(unique(strsplit(a2, "-", fixed = TRUE)[[1]]), "none")
  dom_shared <- length(intersect(s1, s2))
  dom_jac <- if (length(union(s1, s2)) == 0) NA_real_ else dom_shared / length(union(s1, s2))
  dom_class <- if (is.na(dom_match)) NA_character_ else
               if (dom_match) "match" else
               if (dom_shared > 0) "partial_mismatch" else "disjoint_mismatch"

  ## --- risk flags ----------------------------------------------------------
  fl <- character(0)
  if (has_par || n_ip >= 1) fl <- c(fl, "paralogy")
  if (ct_n == 0) fl <- c(fl, "no_ct_hit")
  if (!is.na(ct_qc) && ct_qc < 0.5) fl <- c(fl, "ct_low_qcov")
  if (!is.na(ct_mg) && ct_mg < 1.5) fl <- c(fl, "ct_ambiguous_best")
  if (!is.na(ct_rec) && !ct_rec)    fl <- c(fl, "ct_non_reciprocal")
  if (cb_n == 0) fl <- c(fl, "no_calb_hit")
  if (!is.na(cb_qc) && cb_qc < 0.5) fl <- c(fl, "calb_low_qcov")
  if (!is.na(cb_mg) && cb_mg < 1.5) fl <- c(fl, "calb_ambiguous_best")
  if (!is.na(cb_rec) && !cb_rec)    fl <- c(fl, "calb_non_reciprocal")
  if (!is.na(dom_class) && dom_class == "disjoint_mismatch") fl <- c(fl, "domain_mismatch_disjoint")
  if (!is.na(dom_class) && dom_class == "partial_mismatch")  fl <- c(fl, "domain_partial_mismatch")
  if (pfam_uninf) fl <- c(fl, "pfam_uninformative")
  if (is.na(cr$ct_locus)) fl <- c(fl, "orthomap_no_ortholog")

  ## --- confidence call ------------------------------------------------------
  dom_bad <- !is.na(dom_class) && dom_class == "disjoint_mismatch"
  both_nohit <- ct_n == 0 && cb_n == 0
  both_nonrec <- (!is.na(ct_rec) && !ct_rec) && (!is.na(cb_rec) && !cb_rec)
  one_nonrec  <- (!is.na(ct_rec) && !ct_rec) || (!is.na(cb_rec) && !cb_rec)
  rec_both    <- isTRUE(ct_rec) && isTRUE(cb_rec)
  if (both_nohit || both_nonrec || dom_bad) {
    conf <- "UNRESOLVED"
  } else if (rec_both && length(fl) == 0) {
    conf <- "HIGH"
  } else if (rec_both && length(fl) == 1) {
    conf <- "MEDIUM"
  } else {
    conf <- "LOW"
  }

  ## --- assignment change (defined on the C. tropicalis locus) --------------
  old_loc <- cr$ct_locus
  chg_type <- if (is.na(old_loc) && is.na(ct_loc_new)) "unchanged_no_ortholog"
         else if (is.na(old_loc) && !is.na(ct_loc_new)) "new_ortholog_found"
         else if (!is.na(old_loc) && is.na(ct_loc_new)) "ortholog_lost"
         else if (old_loc == ct_loc_new) "unchanged"
         else "different_locus"
  changed <- !(chg_type %in% c("unchanged", "unchanged_no_ortholog"))
  # a change is only "confident" if the newly assigned C. tropicalis protein is the
  # reciprocal best hit; several ortho_map gaps are filled only by weak one-way hits
  changed_conf <- changed && isTRUE(ct_rec)

  ## --- C. albicans ID-space bookkeeping ------------------------------------
  cb_loc_old <- if (is.na(cr$calb_acc)) NA_character_ else unname(calb_locus_of_up[cr$calb_acc])
  cb_loc_new <- if (is.na(cb_acc)) NA_character_ else unname(calb_locus_of_xp[cb_acc])
  cb_comparable <- !is.na(cb_loc_old) && !is.na(cb_loc_new)

  res[[i]] <- data.frame(
    rank = cr$rank, locus_tag = cr$locus_tag, product = cr$product,
    af_protein = q,
    has_paralog = has_par, n_af_inparalogs = n_ip,
    best_inparalog_acc      = if (n_ip) ip$sseqid[1]   else NA_character_,
    best_inparalog_pident   = if (n_ip) ip$pident[1]   else NA_real_,
    best_inparalog_bitscore = if (n_ip) ip$bitscore[1] else NA_real_,
    ct_locus_orthomap = old_loc, ct_locus_recomputed = ct_loc_new,
    ct_protein_orthomap = cr$ct_protein, ct_protein_recomputed = ct_acc,
    ct_pident = ct_pid, ct_qcov = ct_qc, ct_bitscore = ct_bs,
    ct_margin_vs_second = ct_mg, ct_n_hits = ct_n,
    ct_reciprocal = ct_rec, ct_reciprocal_target = ct_back,
    ct_product_recomputed = if (is.na(ct_acc)) NA_character_ else unname(ct_prod[ct_acc]),
    ct_pident_orthomap = cr$ct_pident,
    ct_qcov_orthomap = if (is.na(cr$ct_qcov)) NA_real_ else as.numeric(cr$ct_qcov) / 100,
    calb_acc_orthomap = cr$calb_acc, calb_acc_recomputed = cb_acc,
    calb_locus_orthomap = cb_loc_old, calb_locus_recomputed = cb_loc_new,
    calb_id_comparable = cb_comparable,
    calb_pident = cb_pid, calb_qcov = cb_qc, calb_bitscore = cb_bs,
    calb_margin_vs_second = cb_mg, calb_n_hits = cb_n,
    calb_reciprocal = cb_rec, calb_reciprocal_target = cb_back,
    calb_product_recomputed = if (is.na(cb_acc)) NA_character_ else unname(calb_prod[cb_acc]),
    af_pfam_architecture = arch_af, ct_pfam_architecture = arch_ct,
    domain_match = dom_match, domain_match_strict = dom_match_strict,
    domain_shared_n = dom_shared, domain_jaccard = dom_jac,
    domain_mismatch_class = dom_class,
    n_risk_flags = length(fl),
    risk_flags = if (length(fl)) paste(fl, collapse = ";") else "none",
    orthology_confidence = conf,
    assignment_changed = changed, assignment_change_type = chg_type,
    assignment_changed_confident = changed_conf,
    confidence_rule = CONF_RULE)
}

out <- do.call(rbind, res)
out <- out[order(out$rank), ]

# ------------------------------------------------------------- 6. write CSVs ---
f1 <- file.path(OUTTAB, paste0("rev_b4_tier1_orthology_confidence", SUFFIX, ".csv"))
write.csv(out, f1, row.names = FALSE, na = "NA")
msg("wrote ", f1, " (", nrow(out), " rows, ", ncol(out), " cols)")

pd <- if (length(para_rows)) do.call(rbind, para_rows) else
  data.frame(locus_tag = character(0), af_protein = character(0),
             inparalog_acc = character(0), pident = numeric(0), qcov = numeric(0),
             bitscore = numeric(0), evalue = numeric(0), inparalog_product = character(0))
f2 <- file.path(OUTTAB, paste0("rev_b4_tier1_paralog_detail", SUFFIX, ".csv"))
write.csv(pd, f2, row.names = FALSE, na = "NA")
msg("wrote ", f2, " (", nrow(pd), " rows)")

# ------------------------------------------------------------- 7. console view --
cat("\n==== Tier-1 orthogroup-aware orthology confidence ====\n")
print(out[, c("rank", "locus_tag", "ct_locus_orthomap", "ct_locus_recomputed",
              "ct_reciprocal", "calb_reciprocal", "domain_match",
              "n_risk_flags", "orthology_confidence", "assignment_changed")],
      row.names = FALSE)
cat("\n-- confidence distribution --\n"); print(table(out$orthology_confidence))
cat("\n-- changed assignments --\n")
print(out[out$assignment_changed, c("locus_tag", "ct_locus_orthomap",
          "ct_locus_recomputed", "assignment_change_type", "orthology_confidence")],
      row.names = FALSE)
cat("\n-- LOW / UNRESOLVED --\n")
print(out[out$orthology_confidence %in% c("LOW", "UNRESOLVED"),
          c("locus_tag", "orthology_confidence", "risk_flags")], row.names = FALSE)
cat("\n-- AFUA_2G01220 full row --\n")
r <- out[out$locus_tag == "AFUA_2G01220", ]
for (k in names(r)) cat(sprintf("  %-26s %s\n", k, as.character(r[[k]])))
msg("done")
