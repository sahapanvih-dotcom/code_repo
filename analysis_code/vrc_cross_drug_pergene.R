#!/usr/bin/env Rscript
# Per-gene weights across all three azoles on the same 55 genomes.
# A variant that shifts all three azoles implicates target/efflux biology; one that shifts only
# voriconazole is more likely a fluke. Same mixed model and GRM as the primary analysis.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
W <- file.path(base,"analysis/vrc_gwas"); out <- file.path(base,"doc/vrc_gene_weighting_2026-08")
GENES <- c("cyp51A","hmg1","cdr1B","srbA","cyp51B","atrR","ptaB","nctA","complexI","cox10")
dat <- read.csv(file.path(out,"tables/model_input_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
az  <- read.csv(file.path(base,"output/table/azole_phenotypes_numeric_fumigatus.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
K <- as.matrix(read.delim(file.path(W,"03_structure/grm_snp.tsv"), row.names=1, check.names=FALSE))
K <- K/mean(diag(K)); K <- K + diag(1e-6, nrow(K))
for (g in GENES) dat[[paste0("car_",g)]] <- as.integer(dat[[paste0("rare_",g)]] > 0)
emma <- function(y, X, K) {
  n <- length(y); p <- ncol(X); S <- diag(n) - X %*% solve(crossprod(X), t(X))
  ev <- eigen(S %*% K %*% S, symmetric=TRUE)
  xi <- pmax(ev$values[seq_len(n-p)],1e-9); U <- ev$vectors[,seq_len(n-p),drop=FALSE]
  eta <- as.vector(crossprod(U,y))
  ll <- function(ld){d<-exp(ld); 0.5*((n-p)*(log((n-p)/(2*pi))-1-log(sum(eta^2/(xi+d))))-sum(log(xi+d)))}
  o <- optimize(ll,c(-12,12),maximum=TRUE); d <- exp(o$maximum)
  sg2 <- sum(eta^2/(xi+d))/(n-p); se2 <- d*sg2
  V <- sg2*K+se2*diag(n); Vi <- solve(V); Ai <- solve(t(X)%*%Vi%*%X)
  b <- Ai%*%t(X)%*%Vi%*%y; list(beta=as.vector(b), se=sqrt(diag(Ai)))
}
res <- lapply(c("voriconazole","itraconazole","posaconazole","amphotericin_B"), function(dr){
  a <- az %>% filter(drug==dr) %>% select(accession, l2=log2_mic)
  m <- dat %>% select(accession, starts_with("car_")) %>% inner_join(a, by="accession")
  m <- m[m$accession %in% rownames(K), ]
  Xg <- as.matrix(m[, paste0("car_", GENES)])
  keep <- apply(Xg, 2, function(v) length(unique(v))>1)
  X <- cbind(1, Xg[, keep, drop=FALSE])
  f <- emma(m$l2, X, K[m$accession, m$accession])
  tibble(drug=dr, gene=GENES[keep], beta=f$beta[-1], se=f$se[-1], n=nrow(m)) %>%
    mutate(p=2*pnorm(-abs(beta/se)), fold=2^beta)
}) %>% bind_rows()
write.csv(res, file.path(out,"tables/cross_drug_per_gene.csv"), row.names=FALSE)
cat("=== per-gene weights across azoles (log2 units)\n")
print(res %>% select(drug, gene, beta, p, fold) %>%
      pivot_wider(names_from=drug, values_from=c(beta,p), names_sep="_") %>% as.data.frame(), digits=3)
cat("\n=== genes with p<0.05 in >1 drug (would indicate consistent biology):\n")
consist <- res %>% filter(p<0.05) %>% count(gene, name="n_drugs") %>% arrange(desc(n_drugs))
print(as.data.frame(consist))
p1 <- ggplot(res %>% filter(drug!="amphotericin_B"),
             aes(beta, reorder(gene, beta), colour=drug)) +
  geom_vline(xintercept=0, linetype=2, colour="grey50") +
  geom_errorbar(aes(xmin=beta-1.96*se, xmax=beta+1.96*se), width=.2, position=position_dodge(.6), linewidth=.45) +
  geom_point(size=2.2, position=position_dodge(.6)) +
  scale_colour_manual(values=c(voriconazole="#B2182B", itraconazole="#33685A", posaconazole="#2166AC")) +
  labs(x="Weight (log2 MIC units - doubling dilutions)", y=NULL, colour=NULL,
       title="Per-gene weights across three azoles, same 55 genomes",
       subtitle="consistent effects across drugs would indicate real target/efflux biology") +
  theme_bw(base_size=11)
ggsave(file.path(out,"figures/cross_drug_per_gene.pdf"), p1, width=9, height=6)
cat("\nwrote cross_drug_per_gene.csv and .pdf\n")
