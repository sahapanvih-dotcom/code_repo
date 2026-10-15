#!/usr/bin/env Rscript
# Summarise the cyp51A positive control: promoter TR genotype + coding substitutions vs voriconazole MIC.
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
m <- read.csv(file.path(base,"output/table/vrc_cyp51A_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
qc <- read.csv(file.path(base,"output/table/vrc_phase0_qc.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
m <- m %>% left_join(qc %>% select(accession, clone_group, group_size), by="accession") %>%
     mutate(TR = factor(promoter_genotype, levels=c("none","TR34","TR46")),
            hotspot = case_when(grepl("L98H",aa_substitutions) ~ "L98H",
                                grepl("Y121F",aa_substitutions) ~ "Y121F",
                                grepl("G54|G138|M220|G448S",aa_substitutions) ~ "other hotspot", TRUE ~ "none"),
            known_mechanism = TR!="none" | hotspot!="none")
cat("=== voriconazole R/S by cyp51A promoter genotype:\n"); print(table(m$TR, m$R_call))
cat("\n=== R without any known cyp51A mechanism:\n")
unex <- m %>% filter(R_call=="R", !known_mechanism) %>% arrange(desc(mic_used))
print(unex %>% select(accession, strain_name, country, mic_used, aa_substitutions, clone_group) %>% as.data.frame())
cat("\n=== S despite a TR genotype (also informative):\n")
print(m %>% filter(R_call=="S", TR!="none") %>% select(accession,strain_name,country,mic_used,TR,aa_substitutions) %>% as.data.frame())
cat("\n=== summary counts:\n")
print(m %>% group_by(known_mechanism) %>% summarise(n=n(), R=sum(R_call=="R"), S=sum(R_call=="S"),
      median_MIC=median(mic_used), .groups="drop") %>% as.data.frame())
ft <- fisher.test(table(m$known_mechanism, m$R_call)); cat("\nFisher exact, known cyp51A mechanism vs VRC R: p =", signif(ft$p.value,3), " OR =", signif(ft$estimate,3), "\n")
write.csv(m, file.path(base,"output/table/vrc_cyp51A_genotypes.csv"), row.names=FALSE, na="")
p <- ggplot(m, aes(TR, log2(mic_used), colour=R_call)) +
  geom_hline(yintercept=0, linetype=2, colour="grey50") +
  geom_jitter(width=.16, height=.06, size=2.6, alpha=.85) +
  scale_colour_manual(values=c(R="#b2182b", S="#2166ac")) +
  scale_y_continuous(breaks=-2:6, labels=2^(-2:6)) +
  labs(x="cyp51A promoter genotype", y="Voriconazole MIC (mg/L)", colour=NULL,
       title="Positive control: cyp51A genotype vs voriconazole MIC",
       subtitle="dashed line = EUCAST breakpoint (1 mg/L); n = 55 A. fumigatus") + theme_bw(base_size=11)
ggsave(file.path(base,"output/figure/vrc_cyp51A_positive_control.pdf"), p, width=7, height=5)
cat("\nwrote output/table/vrc_cyp51A_genotypes.csv and output/figure/vrc_cyp51A_positive_control.pdf\n")
