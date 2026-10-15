#!/usr/bin/env Rscript
# Independence check for every panel gene: does its effect survive among genomes that carry NO known
# cyp51A mechanism, i.e. where the gene would have to act on its own? Any gene whose carriers largely
# overlap cyp51A TR/hotspot genotypes is hitchhiking, not causal (the erg5 lesson, applied a priori).
suppressPackageStartupMessages({library(dplyr)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; out <- file.path(base,"doc/vrc_gene_weighting_2026-08")
GENES <- c("cyp51A","hmg1","cdr1B","srbA","cyp51B","atrR","ptaB","nctA","complexI","cox10")
d <- read.csv(file.path(out,"tables/model_input_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")) %>%
     mutate(km = known_mechanism %in% c("TRUE", TRUE))
res <- lapply(setdiff(GENES,"cyp51A"), function(g){
  car <- d[[paste0("rare_",g)]] > 0
  sub <- d[!d$km, ]; sc <- sub[[paste0("rare_",g)]] > 0
  tb <- table(factor(sc, c(TRUE,FALSE)), factor(sub$R_call, c("R","S")))
  fp <- tryCatch(fisher.test(tb)$p.value, error=function(e) NA)
  lmfit <- tryCatch(coef(summary(lm(log2_mic ~ sc, data=sub))), error=function(e) NULL)
  b  <- if (!is.null(lmfit) && nrow(lmfit)>1) lmfit[2,1] else NA
  bp <- if (!is.null(lmfit) && nrow(lmfit)>1) lmfit[2,4] else NA
  data.frame(gene=g, carriers=sum(car), carriers_with_cyp51A=sum(car & d$km),
             pct_overlap=round(100*sum(car & d$km)/max(sum(car),1)),
             indep_n=nrow(sub), indep_R=tb[1,1], indep_S=tb[1,2],
             indep_fisher_p=round(fp,4), indep_beta_log2=round(b,3), indep_beta_p=round(bp,4),
             verdict=ifelse(is.na(bp) | bp>0.05 | (!is.na(b) && b<0),
                            "not independent of cyp51A", "independent - follow up"))
}) %>% bind_rows() %>% arrange(indep_beta_p)
write.csv(res, file.path(out,"tables/independence_from_cyp51A.csv"), row.names=FALSE, na="")
print(res, row.names=FALSE)
cat("\ngenes with an effect independent of cyp51A:", sum(res$verdict=="independent - follow up"), "of", nrow(res), "\n")
