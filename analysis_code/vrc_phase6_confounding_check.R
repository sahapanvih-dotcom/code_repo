#!/usr/bin/env Rscript
# Confounding check for panel hits: is a candidate variant independent of cyp51A and of lineage?
# A hit whose carriers already carry a cyp51A TR genotype is hitchhiking, not causal. The decisive
# test is to restrict to genomes with NO known cyp51A mechanism, where the candidate must act alone.
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; W <- file.path(base,"analysis/vrc_gwas")
v  <- read.csv(file.path(W,"04_association/panel/panel_variants.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
cy <- read.csv(file.path(base,"output/table/vrc_cyp51A_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")) %>%
      select(accession, strain_name, mic_used, R_call, promoter_genotype, known_mechanism, country)
ps <- read.csv(file.path(base,"output/table/vrc_population_structure.csv"), stringsAsFactors=FALSE) %>% select(accession, clade_k4)
res <- read.csv(file.path(base,"output/table/vrc_candidate_panel_results.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
cand <- res %>% filter(test=="substitution", !is.na(p), p < 0.05, !grepl("^cyp51A", variant)) %>%
        distinct(variant) %>% pull(variant)
base_tab <- cy %>% left_join(ps, by="accession") %>% mutate(known = known_mechanism %in% c("TRUE", TRUE))
out <- lapply(cand, function(vr) {
  g <- sub(":.*","",vr); sb <- sub(".*:","",vr)
  car <- v %>% filter(gene==g) %>% mutate(has = grepl(paste0("(^|;)",sb,"(;|$)"), substitutions)) %>% select(accession, has)
  m <- base_tab %>% left_join(car, by="accession") %>% mutate(has = ifelse(is.na(has), FALSE, has))
  n_car <- sum(m$has); n_car_cyp <- sum(m$has & m$known)
  sub <- m %>% filter(!known)
  tb <- table(factor(sub$has, c(TRUE,FALSE)), factor(sub$R_call, c("R","S")))
  p_ind <- tryCatch(fisher.test(tb)$p.value, error=function(e) NA)
  cl <- table(factor(m$has, c(TRUE,FALSE)), m$clade_k4)
  top_clade <- if (n_car>0) names(which.max(cl[1,])) else NA
  frac_top <- if (n_car>0) round(max(cl[1,])/n_car, 2) else NA
  data.frame(variant=vr, carriers=n_car, carriers_with_cyp51A_mech=n_car_cyp,
             frac_explained_by_cyp51A=round(ifelse(n_car>0, n_car_cyp/n_car, NA),2),
             R_indep=tb[1,1], S_indep=tb[1,2], n_cyp51A_negative=nrow(sub),
             p_independent_of_cyp51A=signif(p_ind,3),
             dominant_clade=top_clade, frac_in_dominant_clade=frac_top,
             verdict=ifelse(is.na(p_ind) | p_ind > 0.05, "confounded / not independent", "independent - follow up"),
             stringsAsFactors=FALSE)
}) %>% bind_rows() %>% arrange(p_independent_of_cyp51A)
write.csv(out, file.path(base,"output/table/vrc_panel_confounding_check.csv"), row.names=FALSE, na="")
print(as.data.frame(out), row.names=FALSE)
cat("\nInterpretation: a candidate whose carriers largely overlap cyp51A TR genotypes, and which loses\n",
    "its association among cyp51A-negative genomes, is hitchhiking on known resistance, not causal.\n")
