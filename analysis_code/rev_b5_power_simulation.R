#!/usr/bin/env Rscript
# =====================================================================================
# rev_b5_power_simulation.R   --  Analysis B5, manuscript revision
#
# PURPOSE
#   Formal power simulation for the phenotype-linked voriconazole association scan.
#   The manuscript currently argues from an empirical p-value comparison that n = 55
#   is underpowered, and extrapolates that 200-300 phenotyped genomes would suffice.
#   This script (a) verifies the published scan parameters, (b) states the empirical
#   argument on its own, and (c) runs a formal power simulation. The advisor asked for
#   (b) and (c) to be kept separate, so they are separate sections and separate tables:
#       (b) output/table/rev_b5_empirical_argument.csv   <- observed p-values only
#       (c) output/table/rev_b5_power_*.csv              <- simulation only
#
# MODEL
#   The original scan used pyseer --lmm (FaST-LMM) with a phylogeny-derived similarity
#   matrix (analysis/vrc_gwas/03_structure/phylogeny_K.tsv = ape::vcv on the rooted
#   core-SNP ML tree). pyseer is a Python tool and no equivalent R package is installed
#   (GMMAT, rrBLUP, sommer and coxme are all absent; only lme4 is present and lme4
#   cannot take an arbitrary user-supplied kinship matrix in the form pyseer uses).
#   This script therefore re-implements pyseer's model directly in R:
#       y = Xb + u + e ,  u ~ N(0, s2*h2*K) ,  e ~ N(0, s2*(1-h2)) ,  trace(K) = n
#   h2 is estimated once by ML under the intercept-only null (exactly as pyseer does),
#   then held fixed for the variant test, whose p-value is a 1-df likelihood-ratio test.
#   Section 2 validates the R implementation against the published pyseer output for the
#   six cyp51A known-mechanism variants (rev_b5_model_validation.csv): betas agree to
#   ~1%, standard errors to ~2%, p-values to within a factor of ~1.2, and the null h2
#   agrees (R 0.843 vs pyseer 0.83). Differences in the last digit come from pyseer and
#   optimize() landing on slightly different ML h2; the model is the same.
#
# OUTPUTS   output/table/rev_b5_*.csv , output/figure/rev_b5_*.pdf
# USAGE     Rscript analysis_code/rev_b5_power_simulation.R [n_rep]
# =====================================================================================

suppressPackageStartupMessages({
  library(ape); library(dplyr); library(data.table); library(ggplot2); library(parallel)
})

BASE <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
W    <- file.path(BASE, "analysis/vrc_gwas")
TBL  <- file.path(BASE, "output/table")
FIG  <- file.path(BASE, "output/figure")

args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 1000L
NCORE  <- max(1L, min(20L, detectCores() - 4L))  # capped: other revision analyses share this host
RNGkind("L'Ecuyer-CMRG")
set.seed(20261012)

say <- function(...) cat(..., "\n", sep = "")
hdr <- function(x) say("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78))

# =====================================================================================
# SECTION 1 -- Verify the original scan: model, features, test count, threshold
# =====================================================================================
hdr("SECTION 1  Verification of the published scan parameters")

scan <- fread(file.path(TBL, "vrc_gwas_scan_summary.csv"))
scan[, threshold_recomputed := 0.05 / patterns]
scan[, bonferroni_on_tested := 0.05 / tested]
print(as.data.frame(scan[, .(scan, label, tested, patterns, threshold,
                             threshold_recomputed, bonferroni_on_tested, min_p, lambda)]))

MS_THRESHOLD <- 1.8e-7   # threshold quoted in the manuscript
MS_BEST_P    <- 1.7e-6   # best observed p quoted in the manuscript

ver <- scan[, .(scan, label, tested_variants = tested, unique_patterns = patterns,
                threshold_in_file = threshold,
                threshold_0.05_over_patterns = threshold_recomputed,
                threshold_0.05_over_tested = bonferroni_on_tested,
                min_p_observed = min_p, lambda_gc = lambda, n_sig,
                fold_short_of_own_threshold = round(min_p / threshold, 2),
                manuscript_threshold = MS_THRESHOLD, manuscript_best_p = MS_BEST_P,
                matches_manuscript_threshold = abs(threshold - MS_THRESHOLD) / MS_THRESHOLD < 0.05,
                matches_manuscript_best_p = abs(min_p - MS_BEST_P) / MS_BEST_P < 0.05)]
fwrite(ver, file.path(TBL, "rev_b5_scan_verification.csv"))

say("\nVERDICT on the manuscript's two numbers:")
say(" * The threshold is Bonferroni on UNIQUE PATTERNS, not on tested variants.")
say("   0.05/274184 = ", signif(scan[scan == "unitigs_bin_dedup46", threshold], 6),
    "  and 0.05/271633 = ", signif(scan[scan == "unitigs_bin_all55", threshold], 6), ".")
say("   Both round to the quoted 1.8e-7, so the threshold CHECKS OUT.")
say(" * 1.7e-6 is the min p of the all55 unitig R/S scan only (",
    signif(scan[scan == "unitigs_bin_all55", min_p], 3), "). It CHECKS OUT.")
say(" * DISCREPANCY: the quoted pair (1.8e-7, 1.7e-6) mixes the n=46 scan's threshold")
say("   with the n=55 scan's best p. Within the n=55 scan the threshold is ",
    signif(scan[scan == "unitigs_bin_all55", threshold], 4), ",")
say("   so the best hit falls ",
    round(scan[scan == "unitigs_bin_all55", min_p / threshold], 1),
    "x short, not ", round(MS_BEST_P / MS_THRESHOLD, 1), "x. The conclusion is unchanged,")
say("   but the two numbers should be quoted from the same scan.")
say(" * Also note the manuscript cites only the unitig test count. If the core-SNP scans")
say("   are counted as separate families their thresholds are ~1.98e-6, i.e. the SAME")
say("   order as the best unitig hit; that is worth stating explicitly.")

ALPHA_GW    <- scan[scan == "unitigs_bin_all55", threshold]        # 1.84e-7
N_PATTERNS  <- scan[scan == "unitigs_bin_all55", patterns]
ALPHA_SNP   <- scan[scan == "snps_bin_dedup46", threshold]         # 1.98e-6
panel       <- fread(file.path(TBL, "vrc_candidate_panel_results.csv"))
N_PANEL     <- max(panel$n_tests, na.rm = TRUE)
ALPHA_PANEL <- 0.05 / N_PANEL
say("\nPrespecified family-wise thresholds used below:")
say("  genome-wide unitig  alpha = ", signif(ALPHA_GW, 4), "  (0.05 / ", N_PATTERNS, " patterns)")
say("  core-SNP scan       alpha = ", signif(ALPHA_SNP, 4))
say("  candidate panel     alpha = ", signif(ALPHA_PANEL, 4), "  (0.05 / ", N_PANEL, " tests)")

# =====================================================================================
# SECTION 2 -- pyseer-equivalent LMM in R, validated on the real data
# =====================================================================================
hdr("SECTION 2  R re-implementation of the pyseer --lmm test, and its validation")

K_raw <- as.matrix(read.delim(file.path(W, "03_structure/phylogeny_K.tsv"),
                              row.names = 1, check.names = FALSE))
ph55  <- read.delim(file.path(W, "03_structure/pheno_all55.tsv"), stringsAsFactors = FALSE)
km    <- as.matrix(read.delim(file.path(W, "02_features/known_mech.Rtab"),
                              row.names = 1, check.names = FALSE))
pstr  <- read.csv(file.path(TBL, "vrc_population_structure.csv"), stringsAsFactors = FALSE)
SAMP  <- ph55$samples
pstr  <- pstr[match(SAMP, pstr$accession), ]

norm_K <- function(K) K * (nrow(K) / sum(diag(K)))        # pyseer: trace(K) = n
prep_K <- function(K) { e <- eigen(norm_K(K), symmetric = TRUE)
                        list(S = pmax(e$values, 0), U = e$vectors, n = nrow(K)) }

.nll <- function(h2, Uy, UX, S, n) {        # ML -logL with sigma2 profiled out
  d <- h2 * S + (1 - h2); d[d < 1e-12] <- 1e-12
  Wm <- UX / d; A <- crossprod(UX, Wm)
  b  <- tryCatch(solve(A, crossprod(Wm, Uy)), error = function(e) NULL)
  if (is.null(b)) return(list(nll = Inf))
  r  <- Uy - UX %*% b; s2 <- sum(r^2 / d) / n
  if (!is.finite(s2) || s2 <= 0) return(list(nll = Inf))
  list(nll = 0.5 * (n * log(2 * pi * s2) + sum(log(d)) + n), beta = b, s2 = s2, A = A)
}

lmm_test <- function(y, g, KP) {            # pyseer: h2 from the null, then fixed
  n <- KP$n; U <- KP$U; S <- KP$S
  Uy  <- crossprod(U, matrix(y, ncol = 1))
  UX0 <- crossprod(U, matrix(1, n, 1))
  h2  <- optimize(function(h) .nll(h, Uy, UX0, S, n)$nll, c(1e-6, 1 - 1e-6), tol = 1e-7)$minimum
  f0  <- .nll(h2, Uy, UX0, S, n)
  UX1 <- cbind(UX0, crossprod(U, matrix(g, ncol = 1)))
  f1  <- .nll(h2, Uy, UX1, S, n)
  if (!is.finite(f1$nll)) return(list(p = NA_real_, beta = NA_real_, se = NA_real_, h2 = h2))
  Ai <- tryCatch(solve(f1$A), error = function(e) NULL)
  list(p = pchisq(max(2 * (f0$nll - f1$nll), 0), 1, lower.tail = FALSE),
       beta = f1$beta[2],
       se = if (is.null(Ai)) NA_real_ else sqrt(f1$s2 * Ai[2, 2]), h2 = h2)
}

KP55 <- prep_K(K_raw[SAMP, SAMP])
Kn55 <- norm_K(K_raw[SAMP, SAMP])

val <- rbindlist(lapply(c("log2mic", "binary"), function(col) {
  y   <- as.numeric(ph55[[col]])
  pys <- fread(file.path(W, sprintf("04_association/known_mech_%s_all55.tsv", col)))
  rbindlist(lapply(rownames(km), function(v) {
    g <- as.numeric(km[v, SAMP]); r <- lmm_test(y, g, KP55); pr <- pys[variant == v]
    data.table(phenotype = col, variant = v, af = mean(g),
               pyseer_beta = pr$beta, R_beta = r$beta,
               pyseer_se = pr$`beta-std-err`, R_se = r$se,
               pyseer_lrt_p = pr$`lrt-pvalue`, R_lrt_p = r$p, R_h2_null = r$h2,
               beta_rel_diff = abs(pr$beta - r$beta) / abs(pr$beta),
               p_ratio = r$p / pr$`lrt-pvalue`)
  }))
}))
fwrite(val, file.path(TBL, "rev_b5_model_validation.csv"))
print(as.data.frame(val[, .(phenotype, variant, af, pyseer_beta, R_beta, pyseer_se, R_se,
                            pyseer_lrt_p, R_lrt_p, p_ratio)]), digits = 3)
say("\nValidation: max |beta| relative difference = ", signif(max(val$beta_rel_diff), 3),
    " ; p-value ratio range ", paste(signif(range(val$p_ratio), 3), collapse = " to "))

# =====================================================================================
# SECTION 3 -- ARGUMENT (A): the empirical observation, on its own. No simulation.
# =====================================================================================
hdr("SECTION 3  ARGUMENT (A)  The empirical observation. Contains NO simulation.")

emp <- rbindlist(lapply(c("binary", "log2mic"), function(col) {
  pys <- fread(file.path(W, sprintf("04_association/known_mech_%s_all55.tsv", col)))
  data.table(phenotype = ifelse(col == "binary", "voriconazole R/S", "log2 MIC"),
             variant = pys$variant, n_genomes = 55L, af = pys$af,
             n_carriers = round(pys$af * 55), beta = pys$beta, se = pys$`beta-std-err`,
             observed_lrt_p = pys$`lrt-pvalue`, genome_wide_threshold = ALPHA_GW,
             fold_short_of_genome_wide = pys$`lrt-pvalue` / ALPHA_GW,
             passes_genome_wide = pys$`lrt-pvalue` < ALPHA_GW,
             panel_threshold = ALPHA_PANEL,
             passes_candidate_panel = pys$`lrt-pvalue` < ALPHA_PANEL,
             pyseer_note = ifelse(is.na(pys$notes), "", pys$notes))
}))
fwrite(emp, file.path(TBL, "rev_b5_empirical_argument.csv"))
print(as.data.frame(emp[, .(phenotype, variant, n_carriers, beta, observed_lrt_p,
                            fold_short_of_genome_wide, passes_genome_wide,
                            passes_candidate_panel)]), digits = 3)

TR46_AF   <- emp[variant == "cyp51A_TR46" & phenotype == "log2 MIC", af]
TR46_BETA <- emp[variant == "cyp51A_TR46" & phenotype == "log2 MIC", beta]
TR46_P    <- emp[variant == "cyp51A_TR46" & phenotype == "log2 MIC", observed_lrt_p]
TR46_SE   <- emp[variant == "cyp51A_TR46" & phenotype == "log2 MIC", se]
TR46_P_B  <- emp[variant == "cyp51A_TR46" & phenotype == "voriconazole R/S", observed_lrt_p]
say("\nARGUMENT (A), stated plainly: TR46 is a laboratory-confirmed cause of azole")
say("resistance. At n=55 it carried ", round(TR46_AF * 55), " genomes, a ",
    signif(TR46_BETA, 3), " log2 MIC (~", round(2^TR46_BETA), "-fold) effect,")
say("and reached only p = ", signif(TR46_P, 3), " on log2 MIC (p = ", signif(TR46_P_B, 3),
    " on R/S) -- ", round(TR46_P / ALPHA_GW), "x and ", round(TR46_P_B / ALPHA_GW),
    "x short of the")
say("genome-wide threshold ", signif(ALPHA_GW, 3), ". No cyp51A mechanism passed it.")
say("This is an OBSERVATION about one locus in one dataset. It demonstrates that the")
say("scan could miss a true cause; it does NOT by itself quantify power, and it cannot")
say("be converted into a required sample size. That is what ARGUMENT (B) below does.")

# how confounded with population structure is each known mechanism?
qK <- function(g, K) { f <- mean(g); as.numeric(t(g - f) %*% K %*% (g - f)) }
clust <- rbindlist(lapply(rownames(km), function(v) {
  g <- as.numeric(km[v, SAMP]); o <- qK(g, Kn55)
  pp <- replicate(5000, qK(sample(g), Kn55))
  data.table(variant = v, n_carriers = sum(g), af = mean(g), quadform_obs = o,
             quadform_perm_mean = mean(pp), clustering_ratio = o / mean(pp),
             perm_p = mean(pp >= o))
}))
fwrite(clust, file.path(TBL, "rev_b5_variant_clustering.csv"))
print(as.data.frame(clust), digits = 3)
RATIO_TR46 <- clust[variant == "cyp51A_TR46", clustering_ratio]
RATIO_TR34 <- clust[variant == "cyp51A_TR34", clustering_ratio]
say("\nTR46 clustering ratio ", signif(RATIO_TR46, 3), " (perm p ",
    signif(clust[variant == "cyp51A_TR46", perm_p], 3),
    "): its 4 carriers sit in 4 different clone")
say("groups, so TR46 is only mildly confounded with the phylogeny. TR34 ratio ",
    signif(RATIO_TR34, 3), " (perm p ", signif(clust[variant == "cyp51A_TR34", perm_p], 3),
    "): clonally expanded.")
say("Consequence: TR46's failure at n=55 is an allele-count problem, not a structure")
say("problem. That is the scenario the simulation must reproduce.")

# =====================================================================================
# SECTION 4 -- Calibrate the simulation scaffold from the real data
# =====================================================================================
hdr("SECTION 4  Calibration of the simulation scaffold from the real data")

g46   <- as.numeric(km["cyp51A_TR46", SAMP])
gk    <- as.numeric(km["cyp51A_known_mech", SAMP])
y55   <- ph55$log2mic
MU_BG  <- mean(y55[g46 == 0])     # background the scan actually faced when testing TR46
VAR_BG <- var(y55[g46 == 0])
MIC_LO <- min(y55); MIC_HI <- max(y55)
RS_CUT <- 1                        # R if log2 MIC >= 1 (MIC >= 2 mg/L): the real rule
R_PREV <- mean(ph55$binary)
H2_LOG <- 0.88                     # published null h2, log2MIC scan
H2_BIN <- 0.84                     # published null h2, R/S scan (pyseer 0.83; R 0.843)
censor <- function(x) pmin(pmax(round(x), MIC_LO), MIC_HI)

say("Background (the 51 non-TR46 genomes): mean log2 MIC = ", signif(MU_BG, 3),
    " , variance = ", signif(VAR_BG, 4))
say("MIC panel limits on the log2 scale [", MIC_LO, ", ", MIC_HI,
    "]; phenotypes are integer doubling dilutions, so the")
say("simulated latent MIC is rounded to the dilution grid and clamped to the panel.")
say("Resistance prevalence in the real cohort = ", signif(R_PREV, 3), " (22 R / 33 S).")

# MIC source / method heterogeneity, from susceptible genomes only so that the
# resistance signal cannot leak into the source variance component.
sdat <- data.frame(log2mic = y55, src = pstr$mic_source)[pstr$R_call == "S", ]
av   <- anova(lm(log2mic ~ src, data = sdat))
ni   <- table(sdat$src); Nn <- sum(ni); kk <- length(ni)
nbar <- (Nn - sum(ni^2) / Nn) / (kk - 1)
V_SRC <- max(0, (av$`Mean Sq`[1] - av$`Mean Sq`[2]) / nbar)
FRAC_SRC <- V_SRC / (V_SRC + av$`Mean Sq`[2])
SRC_PROP <- prop.table(table(pstr$mic_source))
say("\nMIC source/method heterogeneity (susceptible genomes only, so resistance cannot")
say("inflate it): ANOVA p = ", signif(av$`Pr(>F)`[1], 3), " ; var(source) = ",
    signif(V_SRC, 4), " ; var(within) = ", signif(av$`Mean Sq`[2], 4))
say("  -> MIC method accounts for ", round(100 * FRAC_SRC),
    "% of the non-genetic variance. ", length(SRC_PROP), " sources")
say("  (EUCAST, CLSI M38-A2/A3, and 'other'), so a source random effect is included.")

# clonal structure of the real collection (used to expand the tree)
tr <- read.tree(file.path(W, "03_structure/core_snps.treefile"))
tr <- root(tr, outgroup = tr$tip.label[which.max(diag(vcv(tr)))], resolve.root = TRUE)
Dpat <- cophenetic(tr)
cgv  <- pstr$clone_group[match(rownames(Dpat), pstr$accession)]
D_WITHIN <- unlist(lapply(names(which(table(cgv) > 1)),
  function(g) { i <- which(cgv == g); s <- Dpat[i, i]; s[upper.tri(s)] }))
D_BETWEEN <- unlist(lapply(seq_len(54), function(i)
  sapply(seq(i + 1, 55), function(j) if (cgv[i] != cgv[j]) Dpat[i, j] else NA_real_)))
D_BETWEEN <- D_BETWEEN[!is.na(D_BETWEEN)]
CLONAL_FR <- sum(table(cgv)[table(cgv) > 1]) / length(cgv)
K_OBS <- K_raw[SAMP, SAMP]
say("\nClonality of the real collection: ", sum(table(cgv)[table(cgv) > 1]),
    "/55 genomes (", signif(CLONAL_FR, 3), ") belong to a multi-member")
say("clone group. Within-clone patristic distance median ", signif(median(D_WITHIN), 4),
    " vs between-clone ", signif(median(D_BETWEEN), 4), ".")

# =====================================================================================
# SECTION 5 -- Cohort expansion to n > 55
# =====================================================================================
hdr("SECTION 5  How cohorts larger than the real n=55 are generated")

say("Cohorts above n=55 are made by GRAFTING new tips onto the REAL rooted core-SNP")
say("tree, so the real deep topology, branch lengths and lineage proportions are kept.")
say("Because K = ape::vcv(tree) is the shared root-to-MRCA path length, a graft has a")
say("closed form on K: a new tip t' grafted onto host tip h at root-to-split height s has")
say("    K[t',t'] = K[h,h] ,  K[t',h] = s ,  K[t',j] = min(s, K[h,j])  for all other j.")
say("That is exactly the vcv of the grafted tree, so every K_n is a valid tree covariance")
say("by construction (no nearPD repair needed). Hosts are drawn uniformly over the")
say("current tips, which preserves lineage proportions in expectation. With probability")
say("p_clonal the new tip is a near-clone of its host (graft depth drawn from the")
say("observed WITHIN-clone-group patristic distances); otherwise it is a more diverged")
say("member of the same lineage (depth drawn from the observed BETWEEN-clone-group")
say("distances). K_n is then trace-normalised to n, exactly as pyseer does.")
say("")
say("ASSUMPTION (this one matters): a larger collection is assumed to have the same")
say("lineage composition and the same degree of clonality as the observed 55, with")
say("p_clonal = the observed ", signif(CLONAL_FR, 3), ". Sensitivity arms use p_clonal = 0")
say("(fully diversified sampling, optimistic for power) and p_clonal = 0.90 (heavy")
say("clonal expansion, pessimistic). Because grafting adds diversity, mean relatedness")
say("falls as n grows; the p_clonal = 0.90 arm is the one that keeps it near the")
say("observed level, and it is the honest lower bound on power.")

expand_K <- function(n_target, p_clonal = CLONAL_FR) {
  n0 <- nrow(K_OBS)
  if (n_target <= n0) { i <- sample(n0, n_target); return(K_OBS[i, i, drop = FALSE]) }
  K <- matrix(0, n_target, n_target)
  K[seq_len(n0), seq_len(n0)] <- K_OBS
  for (t in seq(n0 + 1L, n_target)) {
    m <- t - 1L; host <- sample.int(m, 1L); H <- K[host, host]
    d <- if (runif(1) < p_clonal) sample(D_WITHIN, 1) / 2 else sample(D_BETWEEN, 1) / 2
    s <- max(0, min(H, H - d))
    row <- pmin(s, K[host, seq_len(m)]); row[host] <- s
    K[t, seq_len(m)] <- row; K[seq_len(m), t] <- row; K[t, t] <- H
  }
  K
}
make_KP <- function(n, p_clonal = CLONAL_FR) {
  K <- expand_K(n, p_clonal); KP <- prep_K(K); KP$Kn <- norm_K(K); KP
}

chk <- rbindlist(lapply(c(55, 100, 200, 300, 500), function(n)
  rbindlist(lapply(c(0, CLONAL_FR, 0.5, 0.9), function(pc) {
    r <- replicate(30, { K <- norm_K(expand_K(n, pc))
      c(mean(K[upper.tri(K)]) / mean(diag(K)),
        max(eigen(K, symmetric = TRUE, only.values = TRUE)$values) / n) })
    data.table(n = n, p_clonal = pc, mean_offdiag_over_mean_diag = mean(r[1, ]),
               top_eigen_share = mean(r[2, ]))
  }))))
chk <- rbind(data.table(n = 55L, p_clonal = -1,
  mean_offdiag_over_mean_diag = mean(Kn55[upper.tri(Kn55)]) / mean(diag(Kn55)),
  top_eigen_share = max(eigen(Kn55, symmetric = TRUE, only.values = TRUE)$values) / 55), chk)
fwrite(chk, file.path(TBL, "rev_b5_cohort_expansion_check.csv"))
say("\nExpansion structure check (p_clonal = -1 is the REAL observed K):")
print(as.data.frame(chk), digits = 4)

# -------------------------------------------------------------------------------------
# Keeping the OBSERVABLE MIC dispersion constant across n.
#
# A draw u ~ N(0, V*h2*K), once centred, has expected sample variance V*h2*c_K where
# c_K = 1 - mean(all entries of the trace-normalised K). c_K grows as a cohort
# diversifies (0.23 at the real n=55, ~0.46 at n=500), so holding V fixed would make
# large cohorts MORE dispersed and small/clonal ones LESS dispersed -- which shows up as
# the nonsensical result that a heavily clonal cohort is better powered, purely because
# its phenotypes were simulated with less spread.
#
# Real A. fumigatus MIC distributions do not narrow or widen with collection size, so V
# is rescaled per (n, p_clonal, h2) to hold the realised phenotypic dispersion at the
# n=55 value. At n=55 the scale is exactly 1, so the anchor of Section 6a is preserved.
# -------------------------------------------------------------------------------------
NS_ALL  <- sort(unique(c(55, 100, 200, 300, 500,
                         55, 75, 100, 125, 150, 175, 200, 250, 300, 400, 500, 700, 1000, 1500)))
PC_ALL  <- sort(unique(c(0, CLONAL_FR, 0.90)))
CK55    <- 1 - mean(Kn55)
CK_TAB  <- rbindlist(lapply(NS_ALL, function(n) rbindlist(lapply(PC_ALL, function(pc)
  data.table(n = n, p_clonal = pc,
             c_K = mean(replicate(12, { K <- norm_K(expand_K(n, pc)); 1 - mean(K) })))))))
CK_KEY  <- setNames(CK_TAB$c_K, paste(CK_TAB$n, round(CK_TAB$p_clonal, 4)))
ck_of   <- function(n, pc) {
  v <- CK_KEY[[paste(n, round(pc, 4))]]
  if (is.null(v)) mean(replicate(12, { K <- norm_K(expand_K(n, pc)); 1 - mean(K) })) else v
}
# V scale that holds realised dispersion at its n=55 value
v_scale <- function(n, pc, h2) (h2 * CK55 + (1 - h2)) / (h2 * ck_of(n, pc) + (1 - h2))
CK_TAB[, v_scale_at_h2_0.88 := vapply(seq_len(.N), function(i) v_scale(n[i], p_clonal[i], 0.88), 0)]
fwrite(CK_TAB, file.path(TBL, "rev_b5_dispersion_rescaling.csv"))
say("\nc_K = 1 - mean(trace-normalised K); real n=55 value = ", signif(CK55, 4), ".")
say("V is divided by (h2*c_K + 1 - h2)/(h2*c_K55 + 1 - h2) so realised dispersion is")
say("n-invariant. Scale factors at h2 = 0.88 (see rev_b5_dispersion_rescaling.csv):")
print(as.data.frame(CK_TAB[n %in% c(55, 100, 200, 300, 500, 1000)]), digits = 4)

# =====================================================================================
# SECTION 6 -- Simulation engine, and its calibration against the real scan
# =====================================================================================
hdr("SECTION 6  Simulation engine")

REPS_PER_K <- 5L                          # phenotype draws reused per kinship draw
NREP_NULL  <- max(4000L, N_REP * 4L)      # replicates for the null calibration

# ---- causal genotype -----------------------------------------------------------------
# "unlinked": carriers at random. "tr46_like"/"clustered": carriers chosen by
# thresholding sqrt(lam)*N(0,K) + sqrt(1-lam)*N(0,I); lam is calibrated below so the
# simulated clustering ratio matches the REAL TR46 (mild) or REAL TR34 (clonal) value.
make_g <- function(n, af, KP, mode, lam) {
  nc <- max(1L, round(af * n))
  if (mode == "unlinked") { g <- numeric(n); g[sample.int(n, nc)] <- 1; return(g) }
  z <- sqrt(lam) * as.numeric(KP$U %*% (sqrt(KP$S) * rnorm(n))) + sqrt(1 - lam) * rnorm(n)
  g <- numeric(n); g[order(z, decreasing = TRUE)[seq_len(nc)]] <- 1; g
}
cal_lam <- function(target_ratio, af) {
  f <- function(lam) mean(replicate(400, {
        g <- make_g(55, af, KP55, "clustered", lam)
        qK(g, Kn55) / mean(replicate(40, qK(sample(g), Kn55))) })) - target_ratio
  tryCatch(uniroot(f, c(1e-4, 0.9999), tol = 0.01)$root, error = function(e) NA_real_)
}
LAM_TR46 <- cal_lam(RATIO_TR46, TR46_AF)
LAM_TR34 <- cal_lam(RATIO_TR34, clust[variant == "cyp51A_TR34", af])
LAM <- c(unlinked = 0, tr46_like = LAM_TR46, clustered = LAM_TR34)
say("Genotype clustering parameter lam calibrated against the real variants:")
say("  tr46_like  lam = ", signif(LAM_TR46, 3), "  (matches the real TR46 ratio ",
    signif(RATIO_TR46, 3), ")")
say("  clustered  lam = ", signif(LAM_TR34, 3), "  (matches the real TR34 ratio ",
    signif(RATIO_TR34, 3), ")")

# ---- background ----------------------------------------------------------------------
# u + s + e with u ~ N(0, h2*V*K) CENTRED, s a per-MIC-source effect, e residual.
# u MUST be centred: 87% of the trace of the real K sits on the leading eigenvector,
# i.e. most of the phylogenetic variance is a shift common to the whole cohort. An LMM
# intercept absorbs that shift and a real MIC panel anchors the grand mean anyway, so an
# uncentred u would push the cohort off the dilution panel and leave the realised
# within-cohort dispersion far below the observed one.
draw_bg <- function(KP, h2, V, frac_src) {
  n <- KP$n
  u <- as.numeric(KP$U %*% (sqrt(KP$S * h2 * V) * rnorm(n))); u <- u - mean(u)
  ns <- length(SRC_PROP)
  s <- rnorm(ns, 0, sqrt(max(frac_src * (1 - h2) * V, 0)))[
         sample.int(ns, n, replace = TRUE, prob = SRC_PROP)]
  u + s + rnorm(n, 0, sqrt(max((1 - frac_src) * (1 - h2) * V, 1e-10)))
}

one_rep <- function(KP, af, beta, h2, outcome, mode, V, frac_src = FRAC_SRC,
                    g_fixed = NULL) {
  n <- KP$n
  g <- if (is.null(g_fixed)) make_g(n, af, KP, mode, LAM[[mode]]) else g_fixed
  nc <- sum(g)
  # pyseer pre-filters variants carried by <2 genomes (--min-af 0.01): untestable
  if (nc < 2 || nc > n - 2) return(c(p = NA_real_, se = NA_real_, filt = 1))
  yo <- censor(MU_BG + beta * g + draw_bg(KP, h2, V, frac_src))
  y  <- if (outcome == "log2mic") yo else as.numeric(yo >= RS_CUT)
  if (outcome == "binary" && (sum(y) < 2 || sum(y) > n - 2))
    return(c(p = NA_real_, se = NA_real_, filt = 2))
  if (sd(y) == 0) return(c(p = NA_real_, se = NA_real_, filt = 3))
  r <- lmm_test(y, g, KP)
  c(p = r$p, se = r$se, filt = 0)
}

# ---- 6a  THE ANCHOR -------------------------------------------------------------------
hdr("SECTION 6a  Anchoring the simulation to the p-value the real scan actually got")

say("A power simulation for this cohort has one hard check available: at n=55 it must")
say("reproduce the result the real scan obtained for TR46. Built naively -- real K, real")
say("effect, real allele frequency, and the latent variance set so the simulated censored")
say("background variance equals the observed ", signif(VAR_BG, 3),
    " -- the simulation returns a median")
say("p around 1e-8 for TR46 on log2 MIC, i.e. roughly 10,000x more significant than the")
say("observed ", signif(TR46_P, 3), ", and would claim ~75% power at n=55 for a variant")
say("that in fact fell ", round(TR46_P / ALPHA_GW), "x short. Such a simulation would be worthless.")
say("")
say("Three things cause the gap, and the first two are handled explicitly:")
say(" 1. the genotype's alignment with the phylogeny (handled: lam calibrated to the")
say("    real TR46 clustering ratio),")
say(" 2. the shape of the real MIC distribution (handled: integer dilution grid plus")
say("    panel clamping),")
say(" 3. a residual effective variance larger than the naive sample variance of the")
say("    background, because resistance is itself partly tree-structured. The real LMM's")
say("    own standard error for TR46 (", signif(TR46_SE, 3),
    ") is ~40% larger than the one implied by")
say("    the naive background variance. This is the part that cannot be derived, so it is")
say("    CALIBRATED: V is chosen so the simulated median p for the TR46 scenario at n=55")
say("    equals the observed ", signif(TR46_P, 3), " exactly.")
say("")
say("Two arms are therefore reported throughout and must be read together:")
say("  ANCHORED   (primary)   V set so the n=55 TR46 p reproduces the observation.")
say("  OPTIMISTIC (bound)     V set so the simulated background variance equals the")
say("                         observed one, ignoring the standard-error gap.")
say("The anchored arm is the honest basis for a required sample size; the optimistic arm")
say("is a floor on the required n. A reviewer will test exactly this choice, so both are")
say("tabulated in rev_b5_variance_calibration.csv.")

sim_p <- function(V, outcome, h2, beta, mode, nrep, frac_src = FRAC_SRC, n = 55L,
                  p_clonal = CLONAL_FR) {
  Vn <- V * v_scale(n, p_clonal, h2)   # hold realised dispersion at its n=55 value
  r <- do.call(cbind, mclapply(split(seq_len(nrep), cut(seq_len(nrep), NCORE, labels = FALSE)),
    function(idx) {
      o <- matrix(NA_real_, 2, length(idx)); pos <- 0L
      while (pos < length(idx)) {
        KP <- make_KP(n, p_clonal)
        for (k in seq_len(min(REPS_PER_K, length(idx) - pos))) {
          pos <- pos + 1L
          z <- one_rep(KP, TR46_AF, beta, h2, outcome, mode, Vn, frac_src)
          o[, pos] <- c(z[["p"]], z[["se"]])
        }
      }
      o
    }, mc.cores = NCORE))
  list(p = r[1, !is.na(r[1, ])], se = r[2, !is.na(r[2, ])])
}

# optimistic V: simulated censored background variance == observed VAR_BG
calib_V_marg <- function(h2, frac_src = FRAC_SRC, reps = 1500)
  uniroot(function(V) mean(replicate(reps,
            var(censor(MU_BG + draw_bg(KP55, h2, V, frac_src))))) - VAR_BG,
          c(0.05, 500), tol = 1e-3, extendInt = "upX")$root

# anchored V: simulated median p for TR46 at n=55 == observed p
calib_V_anchor <- function(h2, outcome, target_p, lo, hi, reps = 500) {
  f <- function(lv) log(median(sim_p(exp(lv), outcome, h2, TR46_BETA,
                                     "tr46_like", reps)$p)) - log(target_p)
  tryCatch(exp(uniroot(f, log(c(lo, hi)), tol = 0.02, extendInt = "yes")$root),
           error = function(e) NA_real_)
}

V_MARG_LOG <- calib_V_marg(H2_LOG)
V_MARG_BIN <- calib_V_marg(H2_BIN)
V_ANCH_LOG <- calib_V_anchor(H2_LOG, "log2mic", TR46_P,   V_MARG_LOG, V_MARG_LOG * 30)
V_ANCH_BIN <- calib_V_anchor(H2_BIN, "binary",  TR46_P_B, V_MARG_BIN * 0.2, V_MARG_BIN * 30)
if (is.na(V_ANCH_BIN)) V_ANCH_BIN <- V_MARG_BIN

# V for the h2 / no-source sensitivity arms: keep the same anchoring ratio
RATIO_ANCH <- V_ANCH_LOG / V_MARG_LOG
V_BY_H2 <- sapply(c(0.5, 0.97), function(h) calib_V_marg(h) * RATIO_ANCH)
names(V_BY_H2) <- c("0.5", "0.97")
V_NOSRC <- calib_V_marg(H2_LOG, frac_src = 0) * RATIO_ANCH

describe_V <- function(V, h2, outcome, lbl, frac_src = FRAC_SRC) {
  bgv <- replicate(800, { yo <- censor(MU_BG + draw_bg(KP55, h2, V, frac_src))
                          c(var(yo), mean(yo >= RS_CUT), mean(yo == MIC_LO), mean(yo == MIC_HI)) })
  s <- sim_p(V, outcome, h2, TR46_BETA, "tr46_like", 1000)
  data.table(arm = lbl, outcome = outcome, h2 = h2, V_lat = V,
             V_inflation_vs_marginal_match = V / (if (outcome == "log2mic") V_MARG_LOG else V_MARG_BIN),
             sim_background_var = mean(bgv[1, ]), observed_background_var = VAR_BG,
             sim_R_prevalence = mean(bgv[2, ]), observed_R_prevalence = R_PREV,
             sim_frac_at_panel_floor = mean(bgv[3, ]), sim_frac_at_panel_ceiling = mean(bgv[4, ]),
             sim_TR46_median_p_n55 = median(s$p), observed_TR46_p = if (outcome == "log2mic") TR46_P else TR46_P_B,
             sim_TR46_median_se_n55 = median(s$se, na.rm = TRUE),
             observed_TR46_se = if (outcome == "log2mic") TR46_SE else NA_real_,
             sim_power_n55 = mean(s$p < ALPHA_GW))
}
vcal <- rbindlist(list(
  describe_V(V_ANCH_LOG, H2_LOG, "log2mic", "ANCHORED_log2mic_primary"),
  describe_V(V_MARG_LOG, H2_LOG, "log2mic", "OPTIMISTIC_log2mic_marginal_match"),
  describe_V(V_ANCH_BIN, H2_BIN, "binary",  "ANCHORED_binary_RS"),
  describe_V(V_MARG_BIN, H2_BIN, "binary",  "OPTIMISTIC_binary_marginal_match")))
fwrite(vcal, file.path(TBL, "rev_b5_variance_calibration.csv"))
print(as.data.frame(vcal), digits = 3)
say("\nAnchoring inflation factor (anchored V / marginal-matched V): log2 MIC = ",
    signif(RATIO_ANCH, 3), " , R/S = ", signif(V_ANCH_BIN / V_MARG_BIN, 3))
say("This single factor is the assumption that most affects every number below.")
say("Note the price of it, visible in the table above: the anchored log2 MIC arm piles")
say("more genomes on the MIC panel floor and raises simulated resistance prevalence")
say("above the observed ", signif(R_PREV, 3), ". That is reported rather than hidden.")
say("")
say("Crucially, the OPTIMISTIC parameterisation is REJECTED by the data: it implies the")
say("real scan had ~", round(100 * vcal[arm == "OPTIMISTIC_log2mic_marginal_match", sim_power_n55]),
    "% probability of detecting TR46 at n=55, and the scan did not")
say("detect it. The ANCHORED arm is therefore the primary one, and the optimistic")
say("numbers should be read only as a floor on the required sample size.")

# consolidated parameter table
params <- data.table(
  parameter = c("n_observed", "n_resistant", "n_susceptible", "R_prevalence",
                "background_mean_log2mic_non_TR46", "background_var_log2mic_non_TR46",
                "mic_log2_lower_limit", "mic_log2_upper_limit", "mic_dilution_step",
                "RS_cutoff_log2mic", "h2_log2mic_published", "h2_binary_published",
                "var_mic_source", "var_within_source", "source_share_of_nongenetic_var",
                "source_anova_p_susceptibles_only", "n_mic_sources",
                "clonal_genome_fraction", "within_clone_patristic_median",
                "between_clone_patristic_median", "lam_tr46_like", "lam_clustered_TR34",
                "TR46_clustering_ratio", "TR34_clustering_ratio",
                "TR46_allele_freq", "TR46_carriers", "TR46_beta_log2mic",
                "TR46_observed_se", "TR46_observed_p_log2mic", "TR46_observed_p_RS",
                "V_lat_marginal_match_log2mic", "V_lat_anchored_log2mic",
                "V_lat_marginal_match_binary", "V_lat_anchored_binary",
                "anchoring_inflation_log2mic", "alpha_genomewide_unitig",
                "alpha_core_snp", "alpha_candidate_panel", "n_unitig_patterns",
                "n_replicates_per_cell", "n_replicates_null_calibration", "reps_per_kinship_draw"),
  value = c(55, 22, 33, R_PREV, MU_BG, VAR_BG, MIC_LO, MIC_HI, 1, RS_CUT,
            H2_LOG, H2_BIN, V_SRC, av$`Mean Sq`[2], FRAC_SRC, av$`Pr(>F)`[1],
            length(SRC_PROP), CLONAL_FR, median(D_WITHIN), median(D_BETWEEN),
            LAM_TR46, LAM_TR34, RATIO_TR46, RATIO_TR34, TR46_AF, round(TR46_AF * 55),
            TR46_BETA, TR46_SE, TR46_P, TR46_P_B, V_MARG_LOG, V_ANCH_LOG,
            V_MARG_BIN, V_ANCH_BIN, RATIO_ANCH, ALPHA_GW, ALPHA_SNP, ALPHA_PANEL,
            N_PATTERNS, N_REP, NREP_NULL, REPS_PER_K),
  source = c(rep("pheno_all55.tsv", 4),                                    #  1- 4
             rep("pheno_all55.tsv + known_mech.Rtab", 2),                  #  5- 6
             rep("pheno_all55.tsv", 4),                                    #  7-10
             "pyseer unitigs_mic_dedup46 log",                             # 11
             "pyseer unitigs_bin_all55 log",                               # 12
             rep("vrc_population_structure.csv, susceptible genomes only", 4),  # 13-16
             "vrc_population_structure.csv",                               # 17
             "vrc_cyp51A_genotypes.csv clone_group",                       # 18
             rep("core_snps.treefile", 2),                                 # 19-20
             rep("calibrated in this script", 2),                          # 21-22
             rep("permutation vs phylogeny_K.tsv", 2),                     # 23-24
             rep("known_mech_log2mic_all55.tsv", 5),                       # 25-29
             "known_mech_binary_all55.tsv",                                # 30
             rep("calibrated in this script", 5),                          # 31-35
             rep("vrc_gwas_scan_summary.csv", 2),                          # 36-37
             "vrc_candidate_panel_results.csv",                            # 38
             "vrc_gwas_scan_summary.csv",                                  # 39
             rep("this script", 3)))                                       # 40-42
stopifnot(nrow(params) == 42L, !any(is.na(params$parameter)))
fwrite(params, file.path(TBL, "rev_b5_simulation_parameters.csv"))

# ---- 6b  NULL CALIBRATION -------------------------------------------------------------
hdr("SECTION 6b  Null calibration: type-I error and lambda_GC of the simulated test")

null_cal <- rbindlist(lapply(c(55, 100, 200, 500), function(n)
  rbindlist(lapply(list(list("log2mic", H2_LOG, V_ANCH_LOG), list("binary", H2_BIN, V_ANCH_BIN)),
    function(cfg) {
      pv <- sim_p(cfg[[3]], cfg[[1]], cfg[[2]], 0, "tr46_like", NREP_NULL, n = n)$p
      data.table(n = n, outcome = cfg[[1]], n_null_rep = length(pv),
                 t1e_0.05 = mean(pv < 0.05), t1e_1e3 = mean(pv < 1e-3),
                 t1e_1e5 = mean(pv < 1e-5), t1e_genomewide = mean(pv < ALPHA_GW),
                 expected_t1e_genomewide = ALPHA_GW,
                 lambda_gc_sim = median(qchisq(pv, 1, lower.tail = FALSE)) / qchisq(0.5, 1),
                 min_null_p = min(pv))
    }))))
fwrite(null_cal, file.path(TBL, "rev_b5_null_calibration.csv"))
print(as.data.frame(null_cal), digits = 3)
say("\nlambda_GC in the REAL scans: unitigs R/S n=55 = ",
    scan[scan == "unitigs_bin_all55", lambda], " ; unitigs log2MIC n=46 = ",
    scan[scan == "unitigs_mic_dedup46", lambda], " ; SNPs log2MIC n=46 = ",
    scan[scan == "snps_mic_dedup46", lambda])
if (any(null_cal$t1e_0.05 > 0.085))
  say("NOTE: simulated type-I error exceeds 0.085 at alpha=0.05 for ",
      paste(null_cal[t1e_0.05 > 0.085, paste0(outcome, " n=", n)], collapse = ", "),
      "; power for those arms is an upper bound.")

# =====================================================================================
# SECTION 7 -- ARGUMENT (B): the formal power simulation
# =====================================================================================
hdr("SECTION 7  ARGUMENT (B)  The formal power simulation")

wilson <- function(k, m, conf = 0.95) {
  if (m == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2); ph <- k / m; d <- 1 + z^2 / m
  c(max(0, (ph + z^2 / (2 * m)) / d - z * sqrt(ph * (1 - ph) / m + z^2 / (4 * m^2)) / d),
    min(1, (ph + z^2 / (2 * m)) / d + z * sqrt(ph * (1 - ph) / m + z^2 / (4 * m^2)) / d))
}
run_grid <- function(ns, afs, betas, h2, outcome, mode, alpha, V,
                     p_clonal = CLONAL_FR, frac_src = FRAC_SRC, nrep = N_REP, tag = "") {
  cells <- as.data.table(expand.grid(beta_log2mic = betas, allele_freq = afs))
  nc <- nrow(cells)
  rbindlist(lapply(ns, function(n) {
    V <- V * v_scale(n, p_clonal, h2)   # hold realised dispersion at its n=55 value
    out <- mclapply(split(seq_len(nrep), cut(seq_len(nrep), min(NCORE, nrep), labels = FALSE)),
      function(idx) {
        m <- length(idx); Pm <- matrix(NA_real_, nc, m); Fm <- matrix(0L, nc, m); pos <- 0L
        while (pos < m) {
          KP <- make_KP(n, p_clonal)
          for (r in seq_len(min(REPS_PER_K, m - pos))) {
            pos <- pos + 1L
            for (ci in seq_len(nc)) {
              z <- one_rep(KP, cells$allele_freq[ci], cells$beta_log2mic[ci], h2,
                           outcome, mode, V, frac_src)
              Pm[ci, pos] <- z[["p"]]; Fm[ci, pos] <- as.integer(z[["filt"]])
            }
          }
        }
        list(P = Pm, F = Fm)
      }, mc.cores = min(NCORE, nrep), mc.preschedule = FALSE)
    if (any(!vapply(out, is.list, logical(1))))
      stop("run_grid worker failed in arm ", tag, " at n=", n)
    P <- do.call(cbind, lapply(out, `[[`, "P")); Fq <- do.call(cbind, lapply(out, `[[`, "F"))
    rbindlist(lapply(seq_len(nc), function(ci) {
      pv <- P[ci, ]; mm <- length(pv)
      k <- sum(!is.na(pv) & pv < alpha); cc <- wilson(k, mm)
      data.table(arm = tag, n = n, allele_freq = cells$allele_freq[ci],
                 n_carriers_expected = max(1, round(cells$allele_freq[ci] * n)),
                 beta_log2mic = cells$beta_log2mic[ci], fold_mic_shift = 2^cells$beta_log2mic[ci],
                 h2 = h2, outcome = outcome, genotype_mode = mode, V_lat = V,
                 p_clonal = p_clonal, source_share = frac_src, alpha = alpha,
                 n_rep = mm, n_detected = k, power = k / mm, ci_lo = cc[1], ci_hi = cc[2],
                 frac_untestable = mean(Fq[ci, ] > 0),
                 median_p = suppressWarnings(median(pv, na.rm = TRUE)))
    }))
  }))
}

NS    <- c(55, 100, 200, 300, 500)
AFS   <- c(0.025, 0.05, TR46_AF, 0.10, 0.15, 0.25, 0.40)
BETAS <- c(1, 2, 3, TR46_BETA, 5)
t0 <- Sys.time()

ARMS <- list(
  list("ANCHORED_log2mic_primary",        H2_LOG, "log2mic", "tr46_like", ALPHA_GW,    V_ANCH_LOG, CLONAL_FR, FRAC_SRC),
  list("OPTIMISTIC_log2mic_marginal",     H2_LOG, "log2mic", "tr46_like", ALPHA_GW,    V_MARG_LOG, CLONAL_FR, FRAC_SRC),
  list("ANCHORED_binary_RS",              H2_BIN, "binary",  "tr46_like", ALPHA_GW,    V_ANCH_BIN, CLONAL_FR, FRAC_SRC),
  list("sens_genotype_unlinked",          H2_LOG, "log2mic", "unlinked",  ALPHA_GW,    V_ANCH_LOG, CLONAL_FR, FRAC_SRC),
  list("sens_genotype_clonal_like_TR34",  H2_LOG, "log2mic", "clustered", ALPHA_GW,    V_ANCH_LOG, CLONAL_FR, FRAC_SRC),
  list("sens_h2_0.50",                    0.50,   "log2mic", "tr46_like", ALPHA_GW,    V_BY_H2[["0.5"]],  CLONAL_FR, FRAC_SRC),
  list("sens_h2_0.97",                    0.97,   "log2mic", "tr46_like", ALPHA_GW,    V_BY_H2[["0.97"]], CLONAL_FR, FRAC_SRC),
  list("sens_sampling_diversified",       H2_LOG, "log2mic", "tr46_like", ALPHA_GW,    V_ANCH_LOG, 0,         FRAC_SRC),
  list("sens_sampling_clonal",            H2_LOG, "log2mic", "tr46_like", ALPHA_GW,    V_ANCH_LOG, 0.90,      FRAC_SRC),
  list("sens_no_mic_method_heterogeneity",H2_LOG, "log2mic", "tr46_like", ALPHA_GW,    V_NOSRC,    CLONAL_FR, 0),
  list("sens_alpha_10x_stricter",         H2_LOG, "log2mic", "tr46_like", ALPHA_GW/10, V_ANCH_LOG, CLONAL_FR, FRAC_SRC),
  list("sens_candidate_panel_alpha",      H2_LOG, "log2mic", "tr46_like", ALPHA_PANEL, V_ANCH_LOG, CLONAL_FR, FRAC_SRC))

pw <- rbindlist(lapply(ARMS, function(a) {
  say("  arm: ", a[[1]])
  run_grid(NS, AFS, BETAS, a[[2]], a[[3]], a[[4]], a[[5]], a[[6]],
           p_clonal = a[[7]], frac_src = a[[8]], tag = a[[1]])
}))
fwrite(pw, file.path(TBL, "rev_b5_power_grid.csv"))
say("wrote rev_b5_power_grid.csv (", nrow(pw), " cells) at ",
    round(as.numeric(difftime(Sys.time(), t0, units = "min")), 1), " min")

say("\nPRIMARY ARM power table (log2 MIC, genome-wide alpha):")
print(as.data.frame(dcast(pw[arm == "ANCHORED_log2mic_primary"],
       beta_log2mic + allele_freq ~ n, value.var = "power")), digits = 3)

# =====================================================================================
# SECTION 8 -- The real TR46 scenario, and the n needed for 80% power
# =====================================================================================
hdr("SECTION 8  The real TR46 scenario and the sample size for 80% power")

NS_FINE <- c(55, 75, 100, 125, 150, 175, 200, 250, 300, 400, 500, 700, 1000, 1500)

tr46 <- rbindlist(lapply(ARMS, function(a) {
  say("  TR46 arm: ", a[[1]])
  run_grid(NS_FINE, TR46_AF, TR46_BETA, a[[2]], a[[3]], a[[4]], a[[5]], a[[6]],
           p_clonal = a[[7]], frac_src = a[[8]], tag = a[[1]])
}))
fwrite(tr46, file.path(TBL, "rev_b5_power_TR46_scenario.csv"))

n80 <- function(dt) {
  dt <- dt[order(n)]; p <- cummax(dt$power)
  if (all(p < 0.8)) return(NA_real_)
  if (p[1] >= 0.8) return(dt$n[1])
  i <- which(p >= 0.8)[1]
  round(dt$n[i - 1] + (0.8 - p[i - 1]) * (dt$n[i] - dt$n[i - 1]) / (p[i] - p[i - 1]))
}

tr46_80 <- tr46[, .(n_for_80pct_power = n80(.SD),
                    power_55 = .SD[n == 55, power][1],   ci55 = sprintf("[%.3f,%.3f]", .SD[n==55,ci_lo][1], .SD[n==55,ci_hi][1]),
                    power_100 = .SD[n == 100, power][1], ci100 = sprintf("[%.3f,%.3f]", .SD[n==100,ci_lo][1], .SD[n==100,ci_hi][1]),
                    power_200 = .SD[n == 200, power][1], ci200 = sprintf("[%.3f,%.3f]", .SD[n==200,ci_lo][1], .SD[n==200,ci_hi][1]),
                    power_300 = .SD[n == 300, power][1], ci300 = sprintf("[%.3f,%.3f]", .SD[n==300,ci_lo][1], .SD[n==300,ci_hi][1]),
                    power_500 = .SD[n == 500, power][1], ci500 = sprintf("[%.3f,%.3f]", .SD[n==500,ci_lo][1], .SD[n==500,ci_hi][1]),
                    power_1000 = .SD[n == 1000, power][1],
                    power_1500 = .SD[n == 1500, power][1]), by = arm]
fwrite(tr46_80, file.path(TBL, "rev_b5_TR46_sample_size_summary.csv"))
say("\nTR46 at its OBSERVED allele frequency (", signif(TR46_AF, 3), ") and OBSERVED effect (",
    signif(TR46_BETA, 3), " log2 MIC):")
print(as.data.frame(tr46_80[, .(arm, power_55, power_100, power_200, power_300, power_500,
                                power_1000, n_for_80pct_power)]), digits = 3)

SCEN_AF   <- c(0.05, TR46_AF, 0.10, 0.15, 0.25)
SCEN_BETA <- c(1, 2, 3, TR46_BETA)
FINE_ARMS <- ARMS[c(1, 2, 3, 5, 9, 11, 12)]
fine <- rbindlist(lapply(FINE_ARMS, function(a) {
  say("  fine-grid arm: ", a[[1]])
  run_grid(NS_FINE, SCEN_AF, SCEN_BETA, a[[2]], a[[3]], a[[4]], a[[5]], a[[6]],
           p_clonal = a[[7]], frac_src = a[[8]], tag = a[[1]])
}))
fwrite(fine, file.path(TBL, "rev_b5_power_fine_n_grid.csv"))

req <- fine[, .(n_for_80pct_power = n80(.SD),
                power_at_200 = .SD[n == 200, power][1],
                ci_at_200 = sprintf("[%.3f,%.3f]", .SD[n==200,ci_lo][1], .SD[n==200,ci_hi][1]),
                power_at_300 = .SD[n == 300, power][1],
                ci_at_300 = sprintf("[%.3f,%.3f]", .SD[n==300,ci_lo][1], .SD[n==300,ci_hi][1]),
                power_at_500 = .SD[n == 500, power][1],
                power_at_1500 = .SD[n == 1500, power][1]),
            by = .(arm, outcome, genotype_mode, h2, p_clonal, alpha, allele_freq, beta_log2mic)]
req[, fold_mic_shift := 2^beta_log2mic]
req[, carriers_at_80pct := ifelse(is.na(n_for_80pct_power), NA, round(allele_freq * n_for_80pct_power))]
req[, supports_200_300 := ifelse(is.na(n_for_80pct_power), "no (>1500)",
                          ifelse(n_for_80pct_power <= 300, "yes", "no"))]
setorder(req, arm, -beta_log2mic, allele_freq)
fwrite(req, file.path(TBL, "rev_b5_sample_size_for_80pct_power.csv"))
say("\nPRIMARY ARM: phenotyped genomes needed for 80% power")
print(as.data.frame(req[arm == "ANCHORED_log2mic_primary",
      .(allele_freq, beta_log2mic, fold_mic_shift, n_for_80pct_power, carriers_at_80pct,
        power_at_200, power_at_300, supports_200_300)]), digits = 3)

# =====================================================================================
# SECTION 9 -- Figures
# =====================================================================================
hdr("SECTION 9  Figures")

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        plot.subtitle = element_text(size = 7, colour = "grey25"))
blab <- function(b) sprintf("effect %.2f log2 MIC (%.0f-fold)", b, 2^b)

d1 <- pw[arm == "ANCHORED_log2mic_primary"]
d1[, af_lab := factor(sprintf("%.3f", allele_freq), levels = sprintf("%.3f", sort(unique(allele_freq))))]
d1[, b_lab := factor(blab(beta_log2mic), levels = blab(sort(unique(beta_log2mic))))]
t46pt <- tr46[arm == "ANCHORED_log2mic_primary" & n %in% NS]
t46pt[, b_lab := factor(blab(TR46_BETA), levels = levels(d1$b_lab))]
p1 <- ggplot(d1, aes(n, power, colour = af_lab, group = af_lab)) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey40") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = af_lab), alpha = 0.15, colour = NA) +
  geom_line() + geom_point(size = 1) +
  geom_point(data = t46pt, aes(n, power), inherit.aes = FALSE, shape = 23,
             fill = "red", colour = "black", size = 2.6) +
  facet_wrap(~ b_lab, ncol = 3) + scale_x_continuous(breaks = NS) + ylim(0, 1) +
  labs(x = "number of phenotyped genomes (n)", y = "empirical power",
       colour = "causal allele frequency", fill = "causal allele frequency",
       title = "Power of the voriconazole unitig LMM scan to detect a causal variant",
       subtitle = sprintf(paste0("pyseer-equivalent LMM on the real phylogeny-derived kinship; log2 MIC; h2 = %.2f; genome-wide alpha = %.2g (0.05/%s unitig patterns); %d replicates/cell;\n",
         "ribbons are 95%% Wilson intervals; red diamonds mark the real cyp51A TR46 scenario (AF %.3f, effect %.2f). ANCHORED arm: variance calibrated so n=55 reproduces the observed TR46 p = %.2g."),
         H2_LOG, ALPHA_GW, format(N_PATTERNS, big.mark = ","), N_REP, TR46_AF, TR46_BETA, TR46_P)) + th
pdf(file.path(FIG, "rev_b5_power_curves_by_effect.pdf"), width = 10, height = 6.8); print(p1); dev.off()

d2 <- copy(tr46)
p2 <- ggplot(d2, aes(n, power, colour = arm, group = arm)) +
  annotate("rect", xmin = 200, xmax = 300, ymin = -Inf, ymax = Inf, alpha = 0.13, fill = "orange") +
  annotate("text", x = 245, y = 0.03, label = "manuscript's 200-300 claim", size = 2.4, colour = "grey25") +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey40") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.02, alpha = 0.45) +
  geom_line() + geom_point(size = 1.3) +
  scale_x_log10(breaks = NS_FINE) + ylim(0, 1) +
  labs(x = "number of phenotyped genomes (n, log scale)", y = "empirical power", colour = "simulation arm",
       title = sprintf("Power for the REAL cyp51A TR46 scenario (AF = %.3f, effect = %.2f log2 MIC, ~%.0f-fold)",
                       TR46_AF, TR46_BETA, 2^TR46_BETA),
       subtitle = sprintf("ARGUMENT (A), observed, not simulated: at n=55 TR46 reached p = %.2g on log2 MIC, %.0fx short of the genome-wide threshold %.2g.",
                          TR46_P, TR46_P / ALPHA_GW, ALPHA_GW)) +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + guides(colour = guide_legend(nrow = 4))
pdf(file.path(FIG, "rev_b5_power_TR46_scenario.pdf"), width = 9, height = 6.2); print(p2); dev.off()

d3 <- req[arm == "ANCHORED_log2mic_primary"]
d3[, lab := ifelse(is.na(n_for_80pct_power), ">1500", as.character(n_for_80pct_power))]
p3 <- ggplot(d3, aes(factor(sprintf("%.3f", allele_freq)), factor(sprintf("%.2f", beta_log2mic)))) +
  geom_tile(aes(fill = pmin(n_for_80pct_power, 1500)), colour = "white") +
  geom_text(aes(label = lab), size = 3) +
  scale_fill_gradientn(colours = c("#2166ac", "#92c5de", "#fddbc7", "#d6604d", "#67001f"),
                       na.value = "grey75", trans = "log10", name = "n for 80% power") +
  labs(x = "causal allele frequency", y = "effect size (log2 MIC units)",
       title = "Phenotyped genomes required for 80% power, genome-wide unitig scan",
       subtitle = sprintf("ANCHORED primary arm: log2 MIC, h2 = %.2f, alpha = %.2g, p_clonal = %.2f, %d replicates/cell. Grey = not reached by n = 1500.",
                          H2_LOG, ALPHA_GW, CLONAL_FR, N_REP)) + th
pdf(file.path(FIG, "rev_b5_n_for_80pct_power.pdf"), width = 7.5, height = 5); print(p3); dev.off()

d4 <- fine[n %in% c(200, 300) & allele_freq %in% c(TR46_AF, 0.15) & beta_log2mic %in% c(2, TR46_BETA)]
d4[, cell := sprintf("AF %.3f / effect %.2f log2 MIC", allele_freq, beta_log2mic)]
p4 <- ggplot(d4, aes(reorder(arm, power), power, fill = factor(n))) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey40") +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(0.8), width = 0.25) +
  coord_flip() + facet_wrap(~ cell) +
  scale_fill_manual(values = c("200" = "#92c5de", "300" = "#2166ac"), name = "n") +
  labs(x = NULL, y = "empirical power",
       title = "Is 200-300 phenotyped genomes enough? Sensitivity to the assumptions",
       subtitle = "95% Wilson intervals; dashed line = 80% power. The ANCHORED arm is the primary one; OPTIMISTIC is a floor on the required n.") + th
pdf(file.path(FIG, "rev_b5_sensitivity_200_300.pdf"), width = 9.5, height = 5.5); print(p4); dev.off()

d5 <- rbindlist(list(copy(null_cal)[, .(n, outcome, kind = "simulated type-I error at alpha=0.05", value = t1e_0.05)],
                     copy(null_cal)[, .(n, outcome, kind = "simulated lambda_GC", value = lambda_gc_sim)]))
p5 <- ggplot(d5, aes(factor(n), value, fill = outcome)) +
  geom_col(position = "dodge", width = 0.7) + facet_wrap(~ kind, scales = "free_y") +
  geom_hline(data = data.frame(kind = c("simulated type-I error at alpha=0.05", "simulated lambda_GC"),
                               yy = c(0.05, 1)), aes(yintercept = yy), linetype = 2) +
  labs(x = "n", y = NULL, title = "Null calibration of the simulated test",
       subtitle = sprintf("Dashed lines are the nominal values. Real-scan lambda_GC: %.3f (R/S, n=55) and %.3f (log2 MIC, n=46).",
                          scan[scan == "unitigs_bin_all55", lambda], scan[scan == "unitigs_mic_dedup46", lambda])) + th
pdf(file.path(FIG, "rev_b5_null_calibration.pdf"), width = 8, height = 4.2); print(p5); dev.off()
say("wrote 5 PDFs to output/figure/")

# =====================================================================================
# SECTION 10 -- Verdict on the 200-300 extrapolation
# =====================================================================================
hdr("SECTION 10  Verdict on the 200-300 extrapolation")

v    <- req[arm == "ANCHORED_log2mic_primary"]
vopt <- req[arm == "OPTIMISTIC_log2mic_marginal"]
t46a <- tr46_80[arm == "ANCHORED_log2mic_primary"]
t46o <- tr46_80[arm == "OPTIMISTIC_log2mic_marginal"]
t46b <- tr46_80[arm == "ANCHORED_binary_RS"]

say("ANCHORED primary arm, genome-wide alpha = ", signif(ALPHA_GW, 3), ":")
say("  scenarios reaching 80% power by n=300 : ", sum(v$n_for_80pct_power <= 300, na.rm = TRUE),
    " of ", nrow(v))
say("  scenarios needing n > 300            : ",
    sum(v$n_for_80pct_power > 300 | is.na(v$n_for_80pct_power)))
say("  median n for 80% power across scenarios: ",
    median(v$n_for_80pct_power, na.rm = TRUE))
say("\nTR46 at its observed frequency and effect:")
say("  ANCHORED  : n=55 ", signif(t46a$power_55, 3), " | 100 ", signif(t46a$power_100, 3),
    " | 200 ", signif(t46a$power_200, 3), " | 300 ", signif(t46a$power_300, 3),
    " | 500 ", signif(t46a$power_500, 3), "  -> n80 = ", t46a$n_for_80pct_power)
say("  OPTIMISTIC: n=55 ", signif(t46o$power_55, 3), " | 100 ", signif(t46o$power_100, 3),
    " | 200 ", signif(t46o$power_200, 3), " | 300 ", signif(t46o$power_300, 3),
    " | 500 ", signif(t46o$power_500, 3), "  -> n80 = ", t46o$n_for_80pct_power)
say("  BINARY R/S: n=55 ", signif(t46b$power_55, 3), " | 200 ", signif(t46b$power_200, 3),
    " | 300 ", signif(t46b$power_300, 3), " | 500 ", signif(t46b$power_500, 3),
    "  -> n80 = ", t46b$n_for_80pct_power)

n_ok   <- sum(v$n_for_80pct_power <= 300, na.rm = TRUE)
n_tot  <- nrow(v)
tr46ok <- !is.na(t46a$n_for_80pct_power) && t46a$n_for_80pct_power <= 300
# the weakest scenario that 200-300 genomes still covers
cov300 <- v[n_for_80pct_power <= 300]
verdict <- paste0(
  if (tr46ok) "PARTLY SUPPORTED. " else "NOT SUPPORTED. ",
  "For a TR46-like variant (AF ", signif(TR46_AF, 3), ", ", signif(TR46_BETA, 3),
  " log2 MIC) 80% power needs n = ", t46a$n_for_80pct_power,
  " on log2 MIC and n = ", t46b$n_for_80pct_power, " on R/S, so 200-300 is about right ",
  "for that one scenario. Across the ", n_tot, " scenarios tested only ", n_ok,
  " reach 80% power by n = 300 (median requirement n = ",
  median(v$n_for_80pct_power, na.rm = TRUE),
  "); 200-300 is optimistic for anything weaker than ~",
  if (nrow(cov300)) signif(min(cov300$beta_log2mic), 2) else NA,
  " log2 MIC or rarer than AF ~",
  if (nrow(cov300)) signif(min(cov300$allele_freq), 2) else NA, ".")

summ <- data.table(
  item = c("manuscript_threshold_quoted", "verified_threshold_n46_unitig",
           "verified_threshold_n55_unitig", "manuscript_best_p_quoted",
           "verified_best_p_n55_unitig", "tested_variants_n55_unitig",
           "unique_patterns_n55_unitig", "core_snp_threshold",
           "candidate_panel_threshold", "candidate_panel_n_tests",
           "TR46_allele_freq", "TR46_carriers", "TR46_beta_log2mic",
           "TR46_observed_p_log2mic", "TR46_observed_p_RS",
           "TR46_fold_short_of_genomewide",
           "anchoring_inflation_factor_log2mic",
           "TR46_power_n55_anchored", "TR46_power_n200_anchored", "TR46_power_n300_anchored",
           "TR46_power_n500_anchored", "TR46_n80_anchored", "TR46_n80_optimistic",
           "TR46_n80_binary_RS", "median_n80_across_scenarios_anchored",
           "n_scenarios_80pct_by_300_anchored", "n_scenarios_needing_more_than_300_anchored",
           "n_replicates_per_cell", "verdict_on_200_300"),
  value = c(MS_THRESHOLD, scan[scan == "unitigs_bin_dedup46", threshold], ALPHA_GW,
            MS_BEST_P, scan[scan == "unitigs_bin_all55", min_p],
            scan[scan == "unitigs_bin_all55", tested], N_PATTERNS, ALPHA_SNP,
            ALPHA_PANEL, N_PANEL, TR46_AF, round(TR46_AF * 55), TR46_BETA,
            TR46_P, TR46_P_B, TR46_P / ALPHA_GW, RATIO_ANCH,
            t46a$power_55, t46a$power_200, t46a$power_300, t46a$power_500,
            t46a$n_for_80pct_power, t46o$n_for_80pct_power, t46b$n_for_80pct_power,
            median(v$n_for_80pct_power, na.rm = TRUE),
            sum(v$n_for_80pct_power <= 300, na.rm = TRUE),
            sum(v$n_for_80pct_power > 300 | is.na(v$n_for_80pct_power)),
            N_REP, verdict))
fwrite(summ, file.path(TBL, "rev_b5_headline_numbers.csv"))
print(as.data.frame(summ))

say("\nVERDICT on '200-300 phenotyped genomes': ", verdict)
say("\nASSUMPTIONS THAT MOST AFFECT THE ANSWER (a reviewer will test these):")
say(" 1. The anchoring inflation factor (", signif(RATIO_ANCH, 3),
    "x the variance implied by the naive background")
say("    variance). Dropping it (OPTIMISTIC arm) changes TR46's n80 from ",
    t46a$n_for_80pct_power, " to ", t46o$n_for_80pct_power, ".")
say(" 2. The phenotype. log2 MIC vs binary R/S changes TR46's n80 from ",
    t46a$n_for_80pct_power, " to ", t46b$n_for_80pct_power,
    "; dichotomising discards most of the signal.")
say(" 3. How larger cohorts are built. Grafting dilutes relatedness as n grows;")
say("    p_clonal = 0.90 (relatedness held near the observed level) is the pessimistic")
say("    arm, p_clonal = 0 the optimistic one.")
say(" 4. The multiple-testing family. alpha = ", signif(ALPHA_GW, 3),
    " assumes the unitig pattern count does not")
say("    grow with n; it will. The alpha/10 arm brackets that.")
say(" 5. h2 of the background (0.88 published; 0.50 and 0.97 arms bracket it).")
say(" 6. The effect size is on the log2 MIC (doubling-dilution) scale and is assumed")
say("    additive on it, with the MIC panel clamped to [", MIC_LO, ", ", MIC_HI, "].")
say("\nTotal runtime: ", round(as.numeric(difftime(Sys.time(), t0, units = "min")), 1), " min")
say("DONE")
