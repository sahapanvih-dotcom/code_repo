#!/usr/bin/env Rscript
# Variance partitioning + per-gene weights for voriconazole log2(MIC), 55 A. fumigatus genomes.
#
# Model:  log2(MIC) = Xb + u + e ,  u ~ N(0, sg2*K) ,  e ~ N(0, se2*I)
#   X = intercept + carrier status for each panel gene (+ optional testing-method covariate)
#   K = genomic relationship matrix (GCTA-style, centred and scaled over 207,707 core SNPs).
#       NOTE: the phylogenetic covariance (ape::vcv on the ML tree) was tried first and is UNSUITABLE
#       here - its mean off-diagonal is 0.77 (nearly an all-ones matrix, because every tip shares most
#       of the root path), which makes the random term absorb the mean and returns a degenerate
#       h2 = 0.96. The SNP-based GRM has mean off-diagonal -0.02 and is the correct choice.
# REML fitted by the EMMA profile-likelihood approach (implemented here; rrBLUP/sommer unavailable).
# Weights are GLS estimates in log2 units, i.e. DOUBLING DILUTIONS of MIC.
# Ridge (glmnet, alpha=0) is fitted as a stability check; bootstrap gives sign consistency.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(glmnet)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
W <- file.path(base,"analysis/vrc_gwas"); out <- file.path(base,"doc/vrc_gene_weighting_2026-08")
set.seed(1)
GENES <- c("cyp51A","hmg1","cdr1B","srbA","cyp51B","atrR","ptaB","nctA","complexI","cox10")

dat <- read.csv(file.path(out,"tables/model_input_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
K   <- as.matrix(read.delim(file.path(W,"03_structure/grm_snp.tsv"), row.names=1, check.names=FALSE))
dat <- dat[dat$accession %in% rownames(K), ]; K <- K[dat$accession, dat$accession]
K <- K / mean(diag(K))                                   # scale so diag averages 1
K <- K + diag(1e-6, nrow(K))                             # ridge for numerical stability
y  <- dat$log2_mic
# carrier status from the RARE encoding (primary): >=1 rare non-synonymous variant or indel
for (g in GENES) dat[[paste0("car_",g)]] <- as.integer(dat[[paste0("rare_",g)]] > 0)
dat$method_class <- ifelse(grepl("EUCAST", dat$mic_method), "EUCAST",
                    ifelse(grepl("CLSI", dat$mic_method), "CLSI", "other"))

# ---------- EMMA-style REML ----------
emma <- function(y, X, K) {
  n <- length(y); p <- ncol(X)
  S <- diag(n) - X %*% solve(crossprod(X), t(X))
  ev <- eigen(S %*% K %*% S, symmetric=TRUE)
  xi <- pmax(ev$values[seq_len(n-p)], 1e-9); U <- ev$vectors[, seq_len(n-p), drop=FALSE]
  eta <- as.vector(crossprod(U, y))
  ll <- function(ld){ d <- exp(ld); 0.5*((n-p)*(log((n-p)/(2*pi))-1-log(sum(eta^2/(xi+d)))) - sum(log(xi+d))) }
  o <- optimize(ll, c(-12,12), maximum=TRUE); d <- exp(o$maximum)
  sg2 <- sum(eta^2/(xi+d))/(n-p); se2 <- d*sg2
  V <- sg2*K + se2*diag(n); Vi <- solve(V)
  A <- t(X) %*% Vi %*% X; Ai <- solve(A); b <- Ai %*% t(X) %*% Vi %*% y
  list(sg2=sg2, se2=se2, delta=d, beta=as.vector(b), se=sqrt(diag(Ai)), logLik=o$objective,
       fitted_fixed=as.vector(X %*% b), X=X)
}
partition <- function(fit, K) {
  Xb <- fit$X[, -1, drop=FALSE] %*% fit$beta[-1]
  vf <- as.numeric(var(Xb)); vg <- fit$sg2 * mean(diag(K)); ve <- fit$se2
  tot <- vf + vg + ve
  c(panel_genes=vf/tot, genetic_background=vg/tot, unexplained=ve/tot, total_variance=tot)
}

X0   <- matrix(1, nrow(dat), 1, dimnames=list(NULL,"(Intercept)"))
Xcyp <- cbind(X0, cyp51A_mechanism = as.integer(dat$known_mechanism %in% c("TRUE", TRUE)))
Xall <- cbind(X0, as.matrix(dat[, paste0("car_", GENES)]))
colnames(Xall) <- c("(Intercept)", GENES)
Xmeth <- cbind(Xall, EUCAST=as.integer(dat$method_class=="EUCAST"), CLSI=as.integer(dat$method_class=="CLSI"))

f0 <- emma(y, X0, K); f1 <- emma(y, Xcyp, K); f2 <- emma(y, Xall, K); f3 <- emma(y, Xmeth, K)
vp <- rbind(`null (lineage only)`=partition(f0,K), `cyp51A mechanism only`=partition(f1,K),
            `full 10-gene panel`=partition(f2,K), `panel + testing method`=partition(f3,K))
vp <- as.data.frame(round(vp, 4))
write.csv(cbind(model=rownames(vp), vp), file.path(out,"tables/variance_partition.csv"), row.names=FALSE)
cat("=== VARIANCE PARTITION (proportion of log2 MIC variance)\n"); print(vp)
cat("\nlineage-only heritability (null model): h2 =", round(f0$sg2*mean(diag(K))/(f0$sg2*mean(diag(K))+f0$se2),3), "\n")

# ---------- weights ----------
wt <- tibble(gene=GENES, beta=f2$beta[-1], se=f2$se[-1]) %>%
  mutate(lo=beta-1.96*se, hi=beta+1.96*se, z=beta/se, p=2*pnorm(-abs(z)),
         fold_change=2^beta, carriers=sapply(GENES, function(g) sum(dat[[paste0("car_",g)]]))) %>%
  arrange(desc(beta))
# ridge (no kinship term; stability check only)
Xr <- as.matrix(dat[, paste0("car_", GENES)]); colnames(Xr) <- GENES
cvr <- cv.glmnet(Xr, y, alpha=0, nfolds=10)
ridge <- as.vector(coef(cvr, s="lambda.min"))[-1]
wt$ridge_beta <- ridge[match(wt$gene, GENES)]
# bootstrap sign stability (OLS on resampled strains)
B <- 1000; sgn <- matrix(NA, B, length(GENES), dimnames=list(NULL, GENES))
for (b in seq_len(B)) {
  i <- sample(nrow(dat), replace=TRUE)
  fit <- try(lm(y[i] ~ Xr[i,,drop=FALSE]), silent=TRUE)
  if (!inherits(fit,"try-error")) sgn[b,] <- sign(coef(fit)[-1])
}
wt$boot_sign_consistency <- sapply(GENES, function(g) {
  s <- sgn[,g]; s <- s[!is.na(s)]; if(!length(s)) NA else round(max(mean(s>0), mean(s<0)),3) })[match(wt$gene, GENES)]
wt <- wt %>% mutate(reliability = ifelse(carriers < 5, "UNRELIABLE (<5 carriers)",
                            ifelse(boot_sign_consistency < 0.8, "unstable sign", "usable")))
# sensitivity: clone-deduplicated subset
rep_acc <- dat %>% group_by(clone_group) %>% slice_max(mic_used, n=1, with_ties=FALSE) %>% pull(accession)
idx <- dat$accession %in% rep_acc
f2d <- emma(y[idx], Xall[idx,,drop=FALSE], K[idx,idx])
wt$beta_dedup46 <- f2d$beta[-1][match(wt$gene, GENES)]
write.csv(wt, file.path(out,"tables/gene_weights.csv"), row.names=FALSE)
cat("\n=== dedup46 sensitivity: variance partition\n"); print(round(partition(f2d, K[idx,idx]),4))
cat("\n=== PER-GENE WEIGHTS (log2 units = doubling dilutions of MIC)\n"); print(as.data.frame(wt), digits=3)

# ---------- collinearity ----------
vif <- diag(solve(cor(Xr)))
cm <- cor(Xr)
write.csv(round(cm,3), file.path(out,"tables/genotype_correlation_matrix.csv"))
cat("\n=== variance inflation factors:\n"); print(round(vif,2))
cat("max |correlation| between genes:", round(max(abs(cm[upper.tri(cm)])),3), "\n")

# ---------- cross-validated incremental R2 ----------
cvR2 <- function(Xm) {
  n <- nrow(Xm); f <- sample(rep(1:10, length.out=n)); pr <- numeric(n)
  for (k in 1:10) { tr <- f!=k
    m <- lm(y[tr] ~ Xm[tr,,drop=FALSE]); cf <- coef(m); cf[is.na(cf)] <- 0
    pr[!tr] <- cbind(1, Xm[!tr,,drop=FALSE]) %*% cf }
  1 - sum((y-pr)^2)/sum((y-mean(y))^2)
}
r2 <- replicate(50, { c(cyp51A_only = cvR2(Xcyp[,-1,drop=FALSE]), full_panel = cvR2(Xr)) })
cat("\n=== cross-validated predictive R2 (50 repeats of 10-fold)\n")
cat("  cyp51A mechanism only:", round(mean(r2["cyp51A_only",]),3), "\n  full 10-gene panel   :", round(mean(r2["full_panel",]),3),
    "\n  increment            :", round(mean(r2["full_panel",])-mean(r2["cyp51A_only",]),3), "\n")
write.csv(data.frame(model=c("cyp51A_only","full_panel"), cv_R2=c(mean(r2["cyp51A_only",]), mean(r2["full_panel",])),
          sd=c(sd(r2["cyp51A_only",]), sd(r2["full_panel",]))), file.path(out,"tables/cv_predictive_r2.csv"), row.names=FALSE)

# ---------- figures ----------
pw <- wt %>% mutate(gene=factor(gene, levels=rev(wt$gene)),
                    sig=ifelse(p<0.05, "p < 0.05", "n.s."))
p1 <- ggplot(pw, aes(beta, gene, colour=sig)) +
  geom_vline(xintercept=0, linetype=2, colour="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=.22, linewidth=.55) +
  geom_point(size=2.8) +
  geom_text(aes(label=paste0("n=",carriers)), hjust=-0.35, vjust=-0.9, size=2.7, colour="grey35") +
  scale_colour_manual(values=c("p < 0.05"="#B2182B","n.s."="#5C6B66")) +
  labs(x="Weight (log2 MIC units - doubling dilutions)", y=NULL, colour=NULL,
       title="Per-gene influence on voriconazole MIC",
       subtitle="mixed model with kinship; bars = 95% CI; n = 55 A. fumigatus genomes") +
  theme_bw(base_size=11)
ggsave(file.path(out,"figures/gene_weights_forest.pdf"), p1, width=8.5, height=5.5)

vpl <- vp %>% mutate(model=rownames(vp)) %>% select(-total_variance) %>%
  pivot_longer(-model, names_to="component", values_to="prop") %>%
  mutate(model=factor(model, levels=rownames(vp)),
         component=factor(component, levels=c("panel_genes","genetic_background","unexplained"),
                          labels=c("Panel genes","Genetic background","Unexplained")))
p2 <- ggplot(vpl, aes(prop, model, fill=component)) + geom_col(width=.62) +
  scale_fill_manual(values=c("Panel genes"="#33685A","Genetic background"="#A8BFD4","Unexplained"="#D9D4CB")) +
  scale_x_continuous(labels=scales::percent) +
  labs(x="Proportion of log2 MIC variance", y=NULL, fill=NULL,
       title="What explains voriconazole MIC variation?",
       subtitle="REML variance partition, 55 A. fumigatus genomes") +
  theme_bw(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(out,"figures/variance_partition.pdf"), p2, width=9, height=4.2)
cat("\nwrote tables and figures to", out, "\n")
