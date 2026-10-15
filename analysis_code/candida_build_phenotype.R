#!/usr/bin/env Rscript
# Phenotype table for the C. tropicalis analysis set.
# Same conventions as the Aspergillus work: numeric MICs only, censored ">x" analysed as 2x,
# "<=x" as x/2. CLSI M27M44S breakpoints for C. tropicalis: fluconazole S <=2, SDD 4, R >=8;
# voriconazole S <=0.12, R >=1. Resistant is coded from the R breakpoint.
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
S <- "/tmp/claude-1000/-home-huanglab/25fa2051-9b57-45bc-894f-755b7f717add/scratchpad/candida"
out <- file.path(base,"doc/Candida")
d <- read.csv(file.path(S,"other_species_mic.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")) %>%
     filter(species=="Candida tropicalis")
BP <- c(fluconazole=8, voriconazole=1, itraconazole=NA, posaconazole=NA,
        amphotericin_b=2, caspofungin=1, micafungin=1, anidulafungin=1, flucytosine=32)
ph <- d %>%
  mutate(drug=gsub(" ","_",tolower(drug)),
         mic_reported=suppressWarnings(as.numeric(mic_mg_per_L)),
         op=ifelse(is.na(mic_operator)|mic_operator=="", "=", mic_operator),
         mic_used=case_when(op %in% c(">",">=") ~ mic_reported*2,
                            op=="<=" ~ mic_reported/2,
                            TRUE ~ mic_reported),
         censored = op != "=") %>%
  filter(!is.na(mic_used), mic_used>0) %>%
  group_by(strain_name, drug) %>% slice(1) %>% ungroup() %>%
  mutate(breakpoint=unname(BP[drug]),
         R_call=ifelse(is.na(breakpoint), NA, ifelse(mic_used>=breakpoint,"R","S")),
         log2_mic=log2(mic_used))
write.csv(ph, file.path(out,"tables/ctropicalis_phenotypes_all_drugs.csv"), row.names=FALSE, na="")
cat("=== per-drug summary\n")
print(ph %>% group_by(drug) %>% summarise(n=n(), R=sum(R_call=="R",na.rm=TRUE), S=sum(R_call=="S",na.rm=TRUE),
      censored=sum(censored), median=median(mic_used), sd_log2=round(sd(log2_mic),2), .groups="drop") %>% as.data.frame())
flc <- ph %>% filter(drug=="fluconazole")
write.csv(flc, file.path(out,"tables/ctropicalis_fluconazole_phenotype.csv"), row.names=FALSE, na="")
cat("\n=== FLUCONAZOLE (primary phenotype)\n")
cat("n =", nrow(flc), " R =", sum(flc$R_call=="R"), " S =", sum(flc$R_call=="S"), "\n")
cat("MIC distribution (mg/L):\n"); print(table(flc$mic_used))
cat("\npublished ERG11 genotype available for:", sum(!is.na(flc$resistance_genotype)), "isolates\n")
cat("carrying Y132F (per publication):", sum(grepl("Y132F", flc$resistance_genotype)), "\n")
cat("published ERG11 copy number stated for:", sum(grepl("copy number", flc$resistance_genotype)), "\n")
cn <- flc %>% filter(grepl("copy number", resistance_genotype)) %>%
  mutate(pub_cn=as.numeric(sub(".*copy number ([0-9]+).*","\\1",resistance_genotype)))
cat("published copy-number distribution:\n"); print(table(cn$pub_cn))
write.csv(cn %>% select(strain_name, biosample, sra_run, mic_used, R_call, resistance_genotype, pub_cn),
          file.path(out,"tables/ctropicalis_published_erg11_cn.csv"), row.names=FALSE, na="")
