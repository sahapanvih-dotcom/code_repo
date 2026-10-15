#!/usr/bin/env Rscript
# Team (mechanism-tier) analysis: group the 10 panel genes into functional teams so that genes with
# too few carriers to test individually are pooled into a testable unit.
#   Sterol pathway : cyp51A, cyp51B, hmg1, erg25B      (the drug target and the pathway feeding it)
#   Efflux         : cdr1B                              (pumps drug out)
#   Regulators     : atrR, ptaB, nctA, srbA             (turn the other two up)
#   Respiration    : complexI, cox10                    (reduce drug uptake)
# Carrier = >=1 rare (AF<=10%) non-synonymous variant or indel in ANY team member.
# The decisive model separates cyp51A from its team, so we can ask whether any team adds anything
# BEYOND the mechanism we already know about.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
W <- file.path(base,"analysis/vrc_gwas"); out <- file.path(base,"doc/vrc_gene_weighting_2026-08")
set.seed(1)
TEAM <- list(sterol_pathway = c("cyp51A","cyp51B","hmg1","erg25B"),
             efflux         = c("cdr1B"),
             regulators     = c("atrR","ptaB","nctA","srbA"),
             respiration    = c("complexI","cox10"))
dat <- read.csv(file.path(out,"tables/model_input_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
K <- as.matrix(read.delim(file.path(W,"03_structure/grm_snp.tsv"), row.names=1, check.names=FALSE))
dat <- dat[dat$accession %in% rownames(K), ]; K <- K[dat$accession, dat$accession]
K <- K/mean(diag(K)); K <- K + diag(1e-6, nrow(K))
have <- function(g) if (paste0("rare_",g) %in% names(dat)) dat[[paste0("rare_",g)]] > 0 else rep(FALSE, nrow(dat))

# team carrier status (any member) and dose (how many members hit)
for (tm in names(TEAM)) {
  mem <- TEAM[[tm]]
  M <- sapply(mem, have); if (is.null(dim(M))) M <- matrix(M, ncol=1)
  dat[[paste0("team_",tm)]]      <- as.integer(rowSums(M) > 0)
  dat[[paste0("teamdose_",tm)]]  <- rowSums(M)
}
# regulators / sterol pathway EXCLUDING cyp51A, so the team can be tested beyond the known mechanism
M <- sapply(setdiff(TEAM$sterol_pathway,"cyp51A"), have)
dat$team_sterol_no_cyp51A <- as.integer(rowSums(M) > 0)
dat$cyp51A_mech <- as.integer(dat$known_mechanism %in% c("TRUE", TRUE))

emma <- function(y, X, K) {
  n <- length(y); p <- ncol(X)
  S <- diag(n) - X %*% solve(crossprod(X), t(X))
  ev <- eigen(S %*% K %*% S, symmetric=TRUE)
  xi <- pmax(ev$values[seq_len(n-p)], 1e-9); U <- ev$vectors[, seq_len(n-p), drop=FALSE]
  eta <- as.vector(crossprod(U, y))
  ll <- function(ld){ d <- exp(ld); 0.5*((n-p)*(log((n-p)/(2*pi))-1-log(sum(eta^2/(xi+d)))) - sum(log(xi+d))) }
  o <- optimize(ll, c(-12,12), maximum=TRUE); d <- exp(o$maximum)
  sg2 <- sum(eta^2/(xi+d))/(n-p); se2 <- d*sg2
  V <- sg2*K + se2*diag(n); Vi <- solve(V); A <- t(X)%*%Vi%*%X; Ai <- solve(A)
  b <- Ai %*% t(X) %*% Vi %*% y
  list(sg2=sg2, se2=se2, beta=as.vector(b), se=sqrt(diag(Ai)), X=X)
}
partition <- function(fit,K){ Xb <- fit$X[,-1,drop=FALSE] %*% fit$beta[-1]
  vf <- as.numeric(var(Xb)); vg <- fit$sg2*mean(diag(K)); ve <- fit$se2; t <- vf+vg+ve
  c(teams=vf/t, background=vg/t, unexplained=ve/t) }
tidyfit <- function(fit, labels, model) tibble(model=model, term=labels, beta=fit$beta[-1], se=fit$se[-1]) %>%
  mutate(lo=beta-1.96*se, hi=beta+1.96*se, fold=2^beta, p=2*pnorm(-abs(beta/se)))

y <- dat$log2_mic; X0 <- matrix(1, nrow(dat), 1)
cat("=== team carrier counts (of", nrow(dat), "genomes)\n")
for (tm in names(TEAM)) cat(sprintf("  %-16s %2d carriers   (members: %s)\n", tm,
    sum(dat[[paste0("team_",tm)]]), paste(TEAM[[tm]], collapse=", ")))
cat(sprintf("  %-16s %2d carriers\n", "sterol (no cyp51A)", sum(dat$team_sterol_no_cyp51A)))

# ---- Model A: four teams as-is
XA <- cbind(X0, as.matrix(dat[, paste0("team_", names(TEAM))])); colnames(XA) <- c("int", names(TEAM))
fA <- emma(y, XA, K)
# ---- Model B: cyp51A separated, teams must add something beyond it
XB <- cbind(X0, cyp51A_mechanism=dat$cyp51A_mech,
            sterol_other=dat$team_sterol_no_cyp51A, efflux=dat$team_efflux,
            regulators=dat$team_regulators, respiration=dat$team_respiration)
fB <- emma(y, XB, K)
# ---- Model C: additive dose across all teams
XC <- cbind(X0, total_dose = rowSums(dat[, paste0("teamdose_", names(TEAM))]))
fC <- emma(y, XC, K)
res <- bind_rows(tidyfit(fA, names(TEAM), "A: four teams"),
                 tidyfit(fB, c("cyp51A mechanism","sterol pathway (excl. cyp51A)","efflux","regulators","respiration"), "B: cyp51A separated"),
                 tidyfit(fC, "total variant dose", "C: additive dose"))
write.csv(res, file.path(out,"tables/team_weights.csv"), row.names=FALSE)
cat("\n=== TEAM WEIGHTS (log2 units = doubling dilutions)\n"); print(as.data.frame(res), digits=3, row.names=FALSE)

vp <- rbind(`A: four teams`=partition(fA,K), `B: cyp51A separated`=partition(fB,K), `C: additive dose`=partition(fC,K))
write.csv(cbind(model=rownames(vp), round(as.data.frame(vp),4)), file.path(out,"tables/team_variance_partition.csv"), row.names=FALSE)
cat("\n=== variance partition\n"); print(round(vp,4))

# ---- independence: teams among cyp51A-negative genomes only
sub <- dat[dat$cyp51A_mech==0, ]
cat("\n=== teams among cyp51A-negative genomes (n =", nrow(sub), ")\n")
ind <- lapply(c(names(TEAM), "sterol_no_cyp51A"), function(tm){
  v <- sub[[paste0("team_",tm)]]
  if (length(unique(v)) < 2) return(data.frame(team=tm, carriers=sum(v), R=NA, S=NA, beta=NA, p=NA, note="no variation"))
  tb <- table(factor(v==1, c(TRUE,FALSE)), factor(sub$R_call, c("R","S")))
  cf <- coef(summary(lm(log2_mic ~ v, data=sub)))
  data.frame(team=tm, carriers=sum(v), R=tb[1,1], S=tb[1,2],
             beta=round(cf[2,1],3), p=round(cf[2,4],4), note="")
}) %>% bind_rows()
write.csv(ind, file.path(out,"tables/team_independence_from_cyp51A.csv"), row.names=FALSE, na="")
print(ind, row.names=FALSE)

# ---- cross-drug robustness (same teams, itraconazole and posaconazole)
az <- read.csv(file.path(base,"output/table/azole_phenotypes_numeric_fumigatus.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
cross <- lapply(c("itraconazole","posaconazole"), function(dr){
  a <- az %>% filter(drug==dr) %>% select(accession, l2=log2_mic)
  m <- dat %>% select(accession, starts_with("team_"), cyp51A_mech) %>% inner_join(a, by="accession")
  Kx <- K[m$accession, m$accession]
  X <- cbind(matrix(1,nrow(m),1), cyp51A_mechanism=m$cyp51A_mech, sterol_other=m$team_sterol_no_cyp51A,
             efflux=m$team_efflux, regulators=m$team_regulators, respiration=m$team_respiration)
  f <- emma(m$l2, X, Kx)
  tidyfit(f, c("cyp51A mechanism","sterol pathway (excl. cyp51A)","efflux","regulators","respiration"), dr) %>% mutate(n=nrow(m))
}) %>% bind_rows()
write.csv(cross, file.path(out,"tables/team_weights_cross_drug.csv"), row.names=FALSE)
cat("\n=== cross-drug check (model B applied to ITC and POS)\n"); print(as.data.frame(cross), digits=3, row.names=FALSE)

# ---- figure
pd <- res %>% filter(model!="C: additive dose") %>%
  mutate(term=factor(term, levels=rev(unique(term))), sig=ifelse(p<0.05,"p < 0.05","n.s."))
p <- ggplot(pd, aes(beta, term, colour=sig)) +
  geom_vline(xintercept=0, linetype=2, colour="grey50") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=.2, linewidth=.55) + geom_point(size=2.8) +
  facet_wrap(~model, ncol=1, scales="free_y") +
  scale_colour_manual(values=c("p < 0.05"="#B2182B","n.s."="#5C6B66")) +
  labs(x="Weight (log2 MIC units - doubling dilutions)", y=NULL, colour=NULL,
       title="Mechanism teams vs voriconazole MIC",
       subtitle="pooling genes with too few carriers to test alone; bars = 95% CI; n = 55") +
  theme_bw(base_size=11)
ggsave(file.path(out,"figures/team_weights.pdf"), p, width=8.5, height=6)
cat("\nwrote team tables and figures\n")
