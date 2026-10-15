#!/usr/bin/env Rscript
# =============================================================================
# rev_b1_essentiality_benchmark.R
# Revision analysis B1 - independent predictive validation.
#
# QUESTION: does within-species evolutionary constraint (piN/piS) predict
# experimental gene essentiality of the C. albicans orthologue BETTER THAN
# simpler baselines (orthologue sequence similarity; gene-level nuisance
# covariates)?
#
# Design, prespecified:
#   outcome   per published C. albicans essentiality SOURCE, separately
#   Model 1   constraint alone               (piN/piS)
#   Model 2   orthologue conservation alone  (% identity / coverage / homolog presence)
#   Model 3   gene-level covariates only     (length, callability, segregating sites,
#                                             annotation availability, orthogroup size)
#   Model 4   Model 1 + Model 2 + Model 3    (integrated)
#   baselines prevalence (constant) and a random predictor
#
# metrics    ROC-AUC and PR-AUC. PR-AUC = AVERAGE PRECISION, computed by
#            summing (recall increment x precision) over the UNIQUE score
#            thresholds, so tied scores are pooled into one operating point
#            (no interpolation, no optimistic tie-breaking). PRROC/precrec are
#            not installed; this is implemented in ap_f() below.
#   CIs      family-blocked paired percentile bootstrap, B = 2000, resampling
#            gene FAMILIES (not genes) with replacement; the same resample
#            index is used for every model, so model deltas are paired.
#   leakage  gene families from single-linkage clustering of the A. fumigatus
#            self-BLAST; grouped (family-blocked) CV; every transform,
#            imputation, scaling and lambda choice is fit inside the training
#            folds only.
#
# HONESTY: every source that yields >= 10 essential and >= 10 non-essential
# eligible genes is reported, whatever the result. Nothing is dropped for being
# unfavourable; skipped sources are reported with the reason.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(glmnet); library(igraph); library(patchwork)
})

base    <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
tabdir  <- file.path(base, "output/table")
figdir  <- file.path(base, "output/figure")
dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

B_BOOT   <- 2000L     # family-blocked bootstrap resamples
K_FOLD   <- 5L        # outer grouped CV folds
N_REPEAT <- 5L        # outer CV repeats
MIN_POS  <- 10L       # minimum essential genes to evaluate a dataset
SEED     <- 20261012L
setDTthreads(0)

say <- function(...) cat(sprintf(...), sep = "")

# ---------------------------------------------------------------- 1. load data
say("== 1. loading inputs\n")
af_raw <- fread(file.path(base, "analysis/target_discovery/constraint_per_gene_dedup.tsv"))
ct_raw <- fread(file.path(base, "analysis/candida/constraint/ctrop_constraint_per_gene.tsv"))
om     <- fread(file.path(base, "analysis/crossgenus/ortho_map.tsv"))
ess    <- fread(file.path(base, "analysis/crossgenus/calb_essentiality.tsv.gz"))
idm    <- fread(file.path(base, "analysis/crossgenus/calb_idmap.tsv"))
selfb  <- fread(file.path(base, "analysis/target_discovery/hom/af293_self.tsv"), header = FALSE)
hum    <- fread(file.path(base, "analysis/target_discovery/hom/af293_vs_human_sens.tsv"), header = FALSE)
# alternative (shrunken) constraint estimator, for the Model-1 sensitivity only
af_alt <- fread(file.path(base, "doc/Aspergillus/target_discovery_2026-08/tables/constraint_all_genes.csv"))[
            , .(locus_tag, piNpiS_shrunk)]
setnames(selfb, c("q", "s", "pident", "alen", "evalue", "bits", "qcov", "scov"))
setnames(hum,   c("q", "s", "pident", "alen", "evalue", "bits", "qcov", "scov"))
say("   AF constraint %d rows | CT constraint %d rows | ortho_map %d | essentiality %d rows\n",
    nrow(af_raw), nrow(ct_raw), nrow(om), nrow(ess))

# ------------------------------------------------- 2. gene families (leakage control)
# Single linkage over NON-SELF self-BLAST hits with pident >= 30%, query coverage
# >= 50% AND subject coverage >= 50%, evalue <= 1e-5. Every A. fumigatus protein is
# a node, so a protein with no qualifying hit is its own (singleton) family.
PID_T <- 30; COV_T <- 50; EV_T <- 1e-5
say("== 2. gene families: single linkage, pident>=%g%%, qcov>=%g%%, scov>=%g%%, E<=%g\n",
    PID_T, COV_T, COV_T, EV_T)
edges <- selfb[q != s & pident >= PID_T & qcov >= COV_T & scov >= COV_T & evalue <= EV_T, .(q, s)]
nodes <- unique(c(af_raw$protein_id, edges$q, edges$s))
g     <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = nodes))
cl    <- components(g)
fam   <- data.table(protein_id = names(cl$membership), family = paste0("FAM", unname(cl$membership)))
fam[, fam_size := .N, by = family]
af_raw <- merge(af_raw, fam, by = "protein_id", all.x = TRUE)
af_raw[is.na(family), `:=`(family = paste0("SING_", protein_id), fam_size = 1L)]
say("   %d edges -> %d families over %d proteins; largest family %d; %d proteins in multi-gene families\n",
    nrow(edges), uniqueN(fam$family), nrow(fam), max(fam$fam_size), sum(fam$fam_size > 1))

# human homolog (sensitive search), best hit per AF protein
humb <- hum[, .(human_pident = max(pident), has_human = 1L), by = q]

# ------------------------------------------------- 3. A. fumigatus predictor table
# Callability / polymorphism-count filter, on this callset's own scale
# (mean_cov_genomes is out of 187 genomes here): >= 150 covered genomes and
# >= 3 segregating sites, and a finite piN/piS.
AF_COV_T <- 150; SEG_T <- 3
af <- copy(af_raw)
af[, tot_seg := syn_seg + nonsyn_seg]
af_elig <- af[mean_cov_genomes >= AF_COV_T & tot_seg >= SEG_T & is.finite(piN_piS)]
say("== 3. AF eligible genes: %d of %d (mean_cov_genomes>=%g, tot_seg>=%d, finite piN/piS)\n",
    nrow(af_elig), nrow(af), AF_COV_T, SEG_T)

unannot_rx <- "uncharacterized|hypothetical|unnamed|unknown function|predicted protein|conserved protein"
afx <- merge(af_elig, af_alt, by = "locus_tag", all.x = TRUE)
afx <- merge(afx, humb, by.x = "protein_id", by.y = "q", all.x = TRUE)
afx[is.na(has_human), `:=`(has_human = 0L, human_pident = 0)]
afx <- merge(afx, om, by.x = "locus_tag", by.y = "af_locus", all.x = TRUE)
AF <- afx[, .(
  gene        = locus_tag,
  family      = family,
  acc         = calb_acc,
  cons        = piN_piS,                                   # PRIMARY constraint estimator
  cons_alt    = piNpiS_shrunk,                             # sensitivity estimator
  calb_pident = calb_pident, calb_qcov = calb_qcov,
  ct_pident   = fifelse(is.na(ct_pident) | ct_protein == "", 0, as.numeric(ct_pident)),
  ct_qcov     = fifelse(is.na(ct_qcov)   | ct_protein == "", 0, as.numeric(ct_qcov)),
  has_ct      = as.integer(!is.na(ct_protein) & ct_protein != ""),
  has_human   = has_human, human_pident = human_pident,
  len         = protein_len,
  callable    = mean_cov_genomes,
  tot_seg     = tot_seg,
  annotated   = as.integer(!grepl(unannot_rx, product, ignore.case = TRUE)),
  og_size     = fam_size)]

# ------------------------------------------------- 4. C. tropicalis predictor table
# CT orthology is A. fumigatus-anchored (ortho_map links AF -> CT and AF -> C. albicans),
# so CT genes inherit the AF gene family as their blocking unit.
CT_CALL_T <- 300
ct <- copy(ct_raw); ct[, tot_seg := syn_seg + nonsyn_seg]
ct_elig <- ct[n_iso_callable >= CT_CALL_T & tot_seg >= SEG_T & is.finite(piNpiS_dedup)]
say("== 4. CT eligible genes: %d of %d (n_iso_callable>=%g, tot_seg>=%d, finite piNpiS_dedup)\n",
    nrow(ct_elig), nrow(ct), CT_CALL_T, SEG_T)
ctx <- merge(om[ct_protein != "", .(af_locus, ct_protein, ct_pident, ct_qcov, calb_acc, calb_pident, calb_qcov)],
             ct_elig, by.x = "ct_protein", by.y = "protein_id")
ctx <- merge(ctx, af_raw[, .(af_locus = locus_tag, family, fam_size)], by = "af_locus", all.x = TRUE)
ctx[is.na(family), `:=`(family = paste0("SING_CT_", ct_protein), fam_size = 1L)]
ctx <- merge(ctx, humb, by.x = "af_locus", by.y = "q", all.x = TRUE)   # AF-anchored human homology
ctx[is.na(has_human), `:=`(has_human = 0L, human_pident = 0)]
CT <- ctx[, .(
  gene        = ct_protein,
  family      = family,
  acc         = calb_acc,
  cons        = piNpiS_dedup,                 # PRIMARY (deduplicated) CT estimator
  cons_alt    = piNpiS_shrunk_dedup,          # sensitivity estimator
  calb_pident = calb_pident, calb_qcov = calb_qcov,
  ct_pident = as.numeric(ct_pident), ct_qcov = as.numeric(ct_qcov), has_ct = 1L,
  has_human = has_human, human_pident = human_pident,
  len         = protein_len,
  callable    = n_iso_callable,
  tot_seg     = tot_seg,
  annotated   = as.integer(!grepl(unannot_rx, product, ignore.case = TRUE)),
  og_size     = fam_size)]

# ------------------------------------------------- 5. outcome definitions
# C. albicans gene IDs -> UniProt accession (the key ortho_map uses).
say("== 5. outcome definitions\n")
setnames(idm, c("acc", "oln", "primary", "syn"))
idl <- idm[, .(id = unlist(strsplit(oln, " "))), by = acc][id != ""]
norm <- function(x) gsub("_", "", sub("^CAALFM_", "", toupper(trimws(x))))
ess[, gid := norm(gene_id)]; idl[, gid := norm(id)]
# GRACEv1_OMeara2015 carries ONLY orf19 names, and the UniProt ordered-locus field
# lists mostly the high-numbered orf19 aliases, so a direct join recovers only 262 of
# its 2,347 genes. A crosswalk is therefore built FROM THIS TABLE ITSELF: the Fu 2021 /
# Fu 2025 rows carry both the assembly-22 systematic ID (gene_id) and the orf19 ID
# (evidence field "orf19=..."), giving verbatim orf19 <-> systematic pairs. No
# identifier is supplied from memory.
cw <- ess[grepl("orf19=orf19", evidence) & grepl("^C[0-9]_", gene_id),
          .(gid_orf19 = norm(sub(".*orf19=(orf19[.0-9]+).*", "\\1", evidence)), gid = norm(gene_id))]
cw <- unique(cw[gid_orf19 != "" & gid != ""])
extra <- unique(merge(cw, idl[, .(acc, gid)], by = "gid")[, .(acc, gid = gid_orf19)])
idl2 <- unique(rbind(idl[, .(acc, gid)], extra))
say("   id crosswalk: %d direct ids + %d orf19 ids recovered from the table itself\n",
    uniqueN(idl$gid), nrow(extra))
e <- merge(ess, idl2, by = "gid", allow.cartesian = TRUE)
say("   essentiality rows mapped to a UniProt accession: %d of %d\n", nrow(e), nrow(ess))

# Assay-condition filters applied BEFORE collapsing to one call per gene:
#  (i) CGD_phenotype: rows from HETEROZYGOUS diploid surveys are dropped. A viable
#      heterozygote is a haploinsufficiency result, not an essentiality result, so
#      those rows are not an admissible NONESSENTIAL observation.
# (ii) CGD_phenotype_indep: additionally drops the CGD rows that re-curate other
#      sources already in this table (PMID 30377286 = Segal 2018,
#      PMID 34764269 = Fu 2021), per README section 3.
# (iii) GRACE "GD" (growth defect) rows arrive as call == UNKNOWN and are always
#      excluded from the binary outcome (counted, not re-binned).
e[, het := source == "CGD_phenotype" & grepl("heterozygous", evidence, ignore.case = TRUE)]
e[, redundant := source == "CGD_phenotype" & grepl("PMID: *(30377286|34764269)", evidence)]

src_defs <- list(
  CGD_phenotype      = quote(source == "CGD_phenotype" & !het),
  CGD_phenotype_indep= quote(source == "CGD_phenotype" & !het & !redundant),
  GRACEv1_OMeara2015 = quote(source == "GRACEv1_OMeara2015"),
  GRACEv1_Fu2021     = quote(source == "GRACEv1_Fu2021"),
  GRACEv2_Fu2021     = quote(source == "GRACEv2_Fu2021"),
  GRACEv3_Fu2025     = quote(source == "GRACEv3_Fu2025"),
  Segal2018          = quote(source == "Segal2018"))

# Two collapse rules, both reported:
#   strict  a gene is kept only if every admissible row in that source agrees
#           (within-source conflicts are EXCLUDED)
#   anyess  a gene is ESSENTIAL if any admissible row says so (conflicts -> essential;
#           this is the rule the earlier manuscript analysis used)
outcomes <- list(); odef <- list()
for (sn in names(src_defs)) {
  sub <- e[eval(src_defs[[sn]])]
  n_unk <- sub[call == "UNKNOWN", .N]
  sub2  <- sub[call %in% c("ESSENTIAL", "NONESSENTIAL")]
  gl <- sub2[, .(nE = sum(call == "ESSENTIAL"), nN = sum(call == "NONESSENTIAL")), by = acc]
  gl[, conflict := nE > 0 & nN > 0]
  strict <- gl[conflict == FALSE, .(acc, y = as.integer(nE > 0))]
  anyess <- gl[, .(acc, y = as.integer(nE > 0))]
  outcomes[[paste(sn, "strict", sep = "|")]] <- strict
  outcomes[[paste(sn, "anyess", sep = "|")]] <- anyess
  odef[[sn]] <- data.table(
    source = sn, rows_in_source = sub[, .N], rows_unknown_dropped = n_unk,
    rows_heterozygous_dropped = if (sn %like% "^CGD") e[source == "CGD_phenotype" & het, .N] else 0L,
    rows_redundant_dropped = if (sn == "CGD_phenotype_indep") e[source == "CGD_phenotype" & redundant & !het, .N] else 0L,
    genes_with_calls = nrow(gl), genes_conflicting = gl[, sum(conflict)],
    genes_strict = nrow(strict), ess_strict = strict[, sum(y)],
    prev_strict = strict[, mean(y)], genes_anyess = nrow(anyess),
    ess_anyess = anyess[, sum(y)], prev_anyess = anyess[, mean(y)])
}
odef <- rbindlist(odef)
print(odef)
fwrite(odef, file.path(tabdir, "rev_b1_outcome_definitions.csv"))

# ------------------------------------------------- 6. metrics
auc_f <- function(s, y) {
  n1 <- sum(y); n0 <- length(y) - n1
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(s)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
# average precision, ties pooled at unique score thresholds
ap_f <- function(s, y) {
  n1 <- sum(y); if (n1 == 0L) return(NA_real_)
  o <- order(s, decreasing = TRUE); ys <- y[o]; ss <- s[o]
  tp <- cumsum(ys); fp <- cumsum(1L - ys)
  keep <- c(ss[-length(ss)] != ss[-1L], TRUE)
  tp <- tp[keep]; fp <- fp[keep]
  P <- tp / (tp + fp); R <- tp / n1
  sum(diff(c(0, R)) * P)
}

# ------------------------------------------------- 7. grouped CV machinery
# feature construction is a pure function of the raw columns (monotone transforms
# only); imputation, centring/scaling and the ridge lambda are learned on the
# training rows of each fold and applied unchanged to the held-out rows.
FSET <- list(
  M1_constraint = c("l_cons"),
  M2_similarity = c("calb_pident", "calb_qcov", "has_ct", "ct_pident", "ct_qcov", "has_human", "human_pident"),
  M3_covariates = c("l_len", "callable", "l_seg", "annotated", "l_og"),
  M4_integrated = c("l_cons", "calb_pident", "calb_qcov", "has_ct", "ct_pident", "ct_qcov",
                    "has_human", "human_pident", "l_len", "callable", "l_seg", "annotated", "l_og"))

make_feats <- function(d, cons_col = "cons") {
  data.table(
    l_cons = -log(pmax(d[[cons_col]], 1e-3) + 1e-3),   # higher = more constrained
    calb_pident = d$calb_pident, calb_qcov = d$calb_qcov,
    has_ct = d$has_ct, ct_pident = d$ct_pident, ct_qcov = d$ct_qcov,
    has_human = d$has_human, human_pident = d$human_pident,
    l_len = log(d$len), callable = d$callable, l_seg = log1p(d$tot_seg),
    annotated = d$annotated, l_og = log(d$og_size))
}

grouped_folds <- function(groups, K, seed) {
  set.seed(seed)
  tb <- table(groups); ug <- sample(names(tb)); sz <- as.integer(tb[ug])
  load <- numeric(K); a <- integer(length(ug))
  ord <- order(-sz)                      # largest families first, greedy balance
  for (i in ord) { k <- which.min(load); a[i] <- k; load[k] <- load[k] + sz[i] }
  unname(setNames(a, ug)[as.character(groups)])
}

oof_predict <- function(X, y, groups, feats, seed) {
  n <- length(y); pred <- rep(NA_real_, n); nrep <- 0L
  acc <- rep(0, n)
  for (rp in seq_len(N_REPEAT)) {
    fo <- grouped_folds(groups, K_FOLD, seed + rp)
    p  <- rep(NA_real_, n); ok <- TRUE
    for (k in seq_len(K_FOLD)) {
      tr <- which(fo != k); te <- which(fo == k)
      if (length(tr) < 20 || length(unique(y[tr])) < 2 || length(te) == 0) { ok <- FALSE; break }
      Xtr <- X[tr, feats, drop = FALSE]; Xte <- X[te, feats, drop = FALSE]
      med <- apply(Xtr, 2, median, na.rm = TRUE); med[!is.finite(med)] <- 0
      for (j in seq_along(feats)) {
        Xtr[is.na(Xtr[, j]), j] <- med[j]; Xte[is.na(Xte[, j]), j] <- med[j]
      }
      mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
      Xtr <- scale(Xtr, mu, sdv); Xte <- scale(Xte, mu, sdv)
      if (length(feats) == 1L) {
        fit <- suppressWarnings(glm.fit(cbind(1, Xtr), y[tr], family = binomial()))
        p[te] <- as.numeric(cbind(1, Xte) %*% fit$coefficients)
      } else {
        fid <- grouped_folds(groups[tr], K_FOLD, seed + 1000L + rp)   # inner CV, same blocking
        cvm <- tryCatch(cv.glmnet(Xtr, y[tr], family = "binomial", alpha = 0,
                                  foldid = fid, type.measure = "deviance"),
                        error = function(z) NULL)
        if (is.null(cvm)) { ok <- FALSE; break }
        p[te] <- as.numeric(predict(cvm, newx = Xte, s = "lambda.min"))
      }
    }
    if (ok && all(is.finite(p))) { acc <- acc + p; nrep <- nrep + 1L }
  }
  if (nrep == 0L) return(rep(NA_real_, n))
  acc / nrep
}

# ------------------------------------------------- 8. run every dataset
# Set REV_B1_FIGURE_ONLY=1 to redraw the figure from the saved result tables without
# refitting (the fit + bootstrap takes ~45 min).
FIG_ONLY <- nzchar(Sys.getenv("REV_B1_FIGURE_ONLY"))
say("== 8. fitting and evaluating%s\n", if (FIG_ONLY) " (SKIPPED: figure-only mode)" else "")
sp_tab <- list(af = AF, ct = CT)
perf <- list(); delt <- list(); skipped <- list()
if (!FIG_ONLY) {

for (ok_ in names(outcomes)) {
  parts <- strsplit(ok_, "|", fixed = TRUE)[[1]]; sn <- parts[1]; rule <- parts[2]
  yt <- outcomes[[ok_]]
  for (sp in names(sp_tab)) {
    d <- merge(sp_tab[[sp]][!is.na(acc) & acc != ""], yt, by = "acc")
    if (nrow(d) == 0 || sum(d$y) < MIN_POS || sum(1 - d$y) < MIN_POS) {
      skipped[[length(skipped) + 1]] <- data.table(source = sn, rule = rule, species = sp,
        n = nrow(d), n_ess = sum(d$y), reason = sprintf("fewer than %d in a class", MIN_POS))
      next
    }
    X  <- as.matrix(make_feats(d))
    Xa <- as.matrix(make_feats(d, "cons_alt"))
    y  <- d$y; gp <- d$family
    scores <- list()
    for (mn in names(FSET)) scores[[mn]] <- oof_predict(X, y, gp, FSET[[mn]], SEED)
    # Model 1 with the alternative (shrunken) estimator - sensitivity only
    if (all(is.finite(Xa[, "l_cons"]))) scores[["M1s_constraint_shrunk"]] <- oof_predict(Xa, y, gp, "l_cons", SEED)
    set.seed(SEED); scores[["B_random"]] <- runif(nrow(d))
    scores[["B_prevalence"]] <- rep(0, nrow(d))
    scores <- scores[!vapply(scores, function(z) all(is.na(z)), TRUE)]
    mods <- names(scores)

    # family-blocked paired bootstrap
    idx_by_fam <- split(seq_len(nrow(d)), gp); nf <- length(idx_by_fam)
    A <- matrix(NA_real_, B_BOOT, length(mods), dimnames = list(NULL, mods))
    Pm <- A
    set.seed(SEED)
    for (b in seq_len(B_BOOT)) {
      ii <- unlist(idx_by_fam[sample.int(nf, nf, replace = TRUE)], use.names = FALSE)
      yb <- y[ii]
      if (sum(yb) < 2L || sum(1L - yb) < 2L) next
      for (m in mods) { sb <- scores[[m]][ii]; A[b, m] <- auc_f(sb, yb); Pm[b, m] <- ap_f(sb, yb) }
    }
    q <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
    for (m in mods) {
      perf[[length(perf) + 1]] <- data.table(
        source = sn, rule = rule, species = sp, model = m,
        n = nrow(d), n_ess = sum(y), prevalence = mean(y),
        n_families = nf,
        roc_auc = auc_f(scores[[m]], y), roc_lo = q(A[, m])[1], roc_hi = q(A[, m])[2],
        pr_auc  = ap_f(scores[[m]], y),  pr_lo  = q(Pm[, m])[1], pr_hi  = q(Pm[, m])[2])
    }
    for (i in seq_along(mods)) for (j in seq_along(mods)) {
      if (i <= j) next
      a <- mods[i]; bb <- mods[j]
      dA <- A[, a] - A[, bb]; dP <- Pm[, a] - Pm[, bb]
      ca <- q(dA); cp <- q(dP)
      delt[[length(delt) + 1]] <- data.table(
        source = sn, rule = rule, species = sp, contrast = paste(a, "-", bb),
        n = nrow(d), prevalence = mean(y),
        d_roc = auc_f(scores[[a]], y) - auc_f(scores[[bb]], y), d_roc_lo = ca[1], d_roc_hi = ca[2],
        d_roc_excl0 = (ca[1] > 0) | (ca[2] < 0),
        d_pr = ap_f(scores[[a]], y) - ap_f(scores[[bb]], y), d_pr_lo = cp[1], d_pr_hi = cp[2],
        d_pr_excl0 = (cp[1] > 0) | (cp[2] < 0),
        boot_p_two_sided_roc = 2 * min(mean(dA <= 0, na.rm = TRUE), mean(dA >= 0, na.rm = TRUE)))
    }
    say("   %-20s %-6s %-2s n=%4d prev=%.3f fam=%4d | AUC M1 %.3f M2 %.3f M3 %.3f M4 %.3f\n",
        sn, rule, sp, nrow(d), mean(y), nf,
        auc_f(scores$M1_constraint, y), auc_f(scores$M2_similarity, y),
        auc_f(scores$M3_covariates, y), auc_f(scores$M4_integrated, y))
  }
}
perf <- rbindlist(perf); delt <- rbindlist(delt)
skipped <- if (length(skipped)) rbindlist(skipped) else data.table()
# every contrast is stored in BOTH directions (the mirror is an exact sign flip of the
# same paired bootstrap, so the CI bounds swap and negate) so that any model pair can
# be read off the table in either orientation
mir <- copy(delt)
mir[, contrast := paste(sub("^(.*) - (.*)$", "\\2", contrast), "-", sub("^(.*) - (.*)$", "\\1", contrast))]
mir[, `:=`(d_roc = -d_roc, d_pr = -d_pr)]
mir[, `:=`(lo2 = -d_roc_hi, hi2 = -d_roc_lo, plo2 = -d_pr_hi, phi2 = -d_pr_lo)]
mir[, `:=`(d_roc_lo = lo2, d_roc_hi = hi2, d_pr_lo = plo2, d_pr_hi = phi2)]
mir[, c("lo2", "hi2", "plo2", "phi2") := NULL]
delt <- rbind(delt, mir)
perf[, species_lab := fifelse(species == "af", "A. fumigatus", "C. tropicalis")]
delt[, species_lab := fifelse(species == "af", "A. fumigatus", "C. tropicalis")]

fwrite(perf, file.path(tabdir, "rev_b1_performance.csv"))
fwrite(delt, file.path(tabdir, "rev_b1_paired_deltas.csv"))
if (nrow(skipped)) fwrite(skipped, file.path(tabdir, "rev_b1_skipped_datasets.csv"))
fam_out <- unique(af_raw[, .(protein_id, locus_tag, family, fam_size)])
fwrite(fam_out[order(-fam_size, family)], file.path(tabdir, "rev_b1_gene_family_blocks.csv"))
} else {
perf <- fread(file.path(tabdir, "rev_b1_performance.csv"))
delt <- fread(file.path(tabdir, "rev_b1_paired_deltas.csv"))
}

say("\n-- primary (strict rule) ROC-AUC --\n")
print(dcast(perf[rule == "strict"], source + species ~ model, value.var = "roc_auc"))
say("\n-- key paired deltas, strict rule --\n")
print(delt[rule == "strict" & contrast %in% c("M4_integrated - M3_covariates",
      "M1_constraint - M2_similarity", "M1_constraint - M3_covariates"),
      .(source, species, contrast, d_roc, d_roc_lo, d_roc_hi, d_roc_excl0)])

# ------------------------------------------------- 9. figure
say("\n== 9. figure\n")
mlab <- c(M1_constraint = "M1 constraint (piN/piS)", M2_similarity = "M2 orthologue similarity",
          M3_covariates = "M3 gene covariates", M4_integrated = "M4 integrated",
          M1s_constraint_shrunk = "M1s constraint (shrunken)",
          B_random = "random", B_prevalence = "prevalence")
mcol <- c("M1 constraint (piN/piS)" = "#1b7a5a", "M2 orthologue similarity" = "#4a7fb5",
          "M3 gene covariates" = "#d98c20", "M4 integrated" = "#6a3d9a",
          "random" = "grey55", "prevalence" = "grey30")
srclab <- c(CGD_phenotype = "CGD curated", CGD_phenotype_indep = "CGD curated (non-redundant)",
            GRACEv1_OMeara2015 = "GRACE v1 O'Meara 2015", GRACEv1_Fu2021 = "GRACE v1 Fu 2021",
            GRACEv2_Fu2021 = "GRACE v2 Fu 2021", GRACEv3_Fu2025 = "GRACE v3 Fu 2025",
            Segal2018 = "Segal 2018 transposon")
th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(colour = "black"), strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8), legend.key.size = unit(3, "mm"),
        legend.position = "bottom", panel.grid.major.y = element_line(colour = "grey92", linewidth = .25))

pp <- perf[rule == "strict" & model %in% c("M1_constraint", "M2_similarity", "M3_covariates",
                                           "M4_integrated", "B_random")]
pp[, model_lab := factor(mlab[model], levels = unname(mlab[c("M1_constraint", "M2_similarity",
                                "M3_covariates", "M4_integrated", "B_random")]))]
pp[, src_lab := factor(srclab[source], levels = unname(srclab))]
pd <- position_dodge(width = .72)

pA <- ggplot(pp, aes(src_lab, roc_auc, colour = model_lab)) +
  geom_hline(yintercept = .5, linetype = "dashed", colour = "grey45", linewidth = .3) +
  geom_errorbar(aes(ymin = roc_lo, ymax = roc_hi), width = .22, position = pd, linewidth = .35) +
  geom_point(size = 1.5, position = pd) +
  facet_wrap(~species_lab, nrow = 1) +
  scale_colour_manual(values = mcol, name = NULL) +
  labs(x = NULL, y = "ROC-AUC (family-blocked 95% CI)",
       title = "A  Discrimination of C. albicans orthologue essentiality") +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

pB <- ggplot(pp, aes(src_lab, pr_auc, colour = model_lab)) +
  geom_point(aes(y = prevalence), shape = 95, size = 5, colour = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = pr_lo, ymax = pr_hi), width = .22, position = pd, linewidth = .35) +
  geom_point(size = 1.5, position = pd) +
  facet_wrap(~species_lab, nrow = 1) +
  scale_colour_manual(values = mcol, name = NULL) +
  labs(x = NULL, y = "PR-AUC / average precision",
       title = "B  PR-AUC (black dash = prevalence baseline)") +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

key <- c("M1_constraint - M2_similarity", "M1_constraint - M3_covariates", "M4_integrated - M3_covariates")
dd <- delt[rule == "strict" & contrast %in% key]
dd[, src_lab := factor(srclab[source], levels = unname(srclab))]
dd[, con_lab := factor(contrast, levels = key,
      labels = c("M1 - M2", "M1 - M3", "M4 - M3"))]
pC <- ggplot(dd, aes(src_lab, d_roc, colour = con_lab)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = .3) +
  geom_errorbar(aes(ymin = d_roc_lo, ymax = d_roc_hi), width = .22, position = pd, linewidth = .35) +
  geom_point(aes(shape = d_roc_excl0), size = 1.6, position = pd) +
  facet_wrap(~species_lab, nrow = 1) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "95% CI excludes 0") +
  scale_colour_manual(values = c("M1 - M2" = "#1b7a5a", "M1 - M3" = "#b2182b", "M4 - M3" = "#6a3d9a"), name = NULL) +
  labs(x = NULL, y = "paired delta ROC-AUC",
       title = "C  Paired model contrasts (same family-blocked resamples)") +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

fig <- pA / pB / pC
ggsave(file.path(figdir, "rev_b1_essentiality_benchmark.pdf"), fig,
       width = 9.2, height = 10.6, device = cairo_pdf)
say("   wrote %s\n", file.path(figdir, "rev_b1_essentiality_benchmark.pdf"))
say("== done\n")
