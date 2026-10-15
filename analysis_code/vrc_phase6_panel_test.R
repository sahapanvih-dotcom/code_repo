#!/usr/bin/env Rscript
# Association test on the candidate-gene panel (few tests -> adequately powered at n=55).
# Two tests per gene: (a) any protein-changing variant vs voriconazole R/S (Fisher),
#                     (b) per-substitution test for substitutions seen in >=3 genomes.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; W <- file.path(base,"analysis/vrc_gwas")
v  <- read.csv(file.path(W,"04_association/panel/panel_variants.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
ph <- read.csv(file.path(base,"output/table/vrc_cyp51A_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")) %>%
      select(accession, strain_name, mic_used, R_call, promoter_genotype, known_mechanism, clone_group, country)
cat("panel call status:\n"); print(table(v$status))
d <- v %>% filter(status=="ok") %>% left_join(ph, by="accession")
# dedup set: one representative per clone group
rep_acc <- ph %>% group_by(clone_group) %>% slice_max(mic_used, n=1, with_ties=FALSE) %>% pull(accession)
test_set <- function(df, label) {
  # (a) gene burden
  g <- df %>% mutate(has = n_subs > 0) %>% group_by(gene) %>%
    summarise(n=n(), R_with=sum(has & R_call=="R"), R_without=sum(!has & R_call=="R"),
              S_with=sum(has & R_call=="S"), S_without=sum(!has & R_call=="S"), .groups="drop") %>%
    rowwise() %>% mutate(p = tryCatch(fisher.test(matrix(c(R_with,R_without,S_with,S_without),2))$p.value, error=function(e) NA)) %>%
    ungroup() %>% mutate(set=label, test="gene_burden", variant=NA) %>% arrange(p)
  # (b) individual substitutions
  subs <- df %>% filter(!is.na(substitutions)) %>%
    separate_rows(substitutions, sep=";") %>% filter(substitutions!="") %>%
    mutate(variant=paste0(gene,":",substitutions)) %>% distinct(accession, variant, R_call, gene)
  cnt <- subs %>% count(variant, name="carriers") %>% filter(carriers>=3, carriers<=nrow(df)/length(unique(df$gene))-3)
  allacc <- unique(df$accession); rs <- df %>% distinct(accession, R_call)
  st <- lapply(cnt$variant, function(vr) {
    car <- subs$accession[subs$variant==vr]
    tb <- table(factor(rs$accession %in% car, c(TRUE,FALSE)), factor(rs$R_call, c("R","S")))
    data.frame(variant=vr, gene=sub(":.*","",vr), n=nrow(rs), R_with=tb[1,1], R_without=tb[2,1], S_with=tb[1,2], S_without=tb[2,2],
               p=fisher.test(tb)$p.value, stringsAsFactors=FALSE)
  })
  st <- if (length(st)) bind_rows(st) %>% mutate(set=label, test="substitution") %>% arrange(p) else NULL
  bind_rows(g, st)
}
res <- bind_rows(test_set(d %>% filter(accession %in% rep_acc), "dedup46"), test_set(d, "all55"))
res <- res %>% group_by(set, test) %>% mutate(n_tests=n(), p_bonf=pmin(1, p*n_tests)) %>% ungroup() %>% arrange(p)
write.csv(res, file.path(base,"output/table/vrc_candidate_panel_results.csv"), row.names=FALSE, na="")
cat("\n=== top results (dedup46):\n")
print(res %>% filter(set=="dedup46") %>% slice_min(p, n=12) %>%
      select(test, gene, variant, R_with, R_without, S_with, S_without, p, n_tests, p_bonf) %>% as.data.frame(), digits=3)
cat("\n=== significant after Bonferroni within test family:\n")
sig <- res %>% filter(p_bonf < 0.05)
print(if (nrow(sig)) as.data.frame(sig %>% select(set,test,gene,variant,R_with,R_without,S_with,S_without,p,p_bonf)) else "none", digits=3)
