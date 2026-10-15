#!/usr/bin/env Rscript
# Combine ERG11 protein genotype + copy number + phenotype into the analysis table.
suppressPackageStartupMessages({library(data.table); library(dplyr)})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; C <- file.path(B,"analysis/candida")
out <- file.path(B,"doc/Candida")
PANEL <- c(ERG11="CTRG_05283", TAC1="CTRG_05307", UPC2="CTRG_03159", MRR1="CTRG_02269",
           MDR1="CTRG_02412", CDR2="CTRG_02337", ERG3="CTRG_04480", ERG5="CTRG_04983",
           ERG6="CTRG_02396", ERG2="CTRG_04378")
have <- sub("\\.vcf\\.gz$","",basename(Sys.glob(file.path(C,"vcf/*.vcf.gz"))))
fs <- file.path(C,"depth", paste0(have,".regions.bed.gz")); fs <- fs[file.exists(fs)]
cn <- rbindlist(lapply(fs, function(f){
  x <- fread(f, col.names=c("c","s","e","g","dp")); med <- median(x$dp[x$dp>0])
  d <- data.table(accession=sub("\\.regions\\.bed\\.gz$","",basename(f)), median_depth=med)
  for (nm in names(PANEL)) d[[paste0("CN_",nm)]] <- round(2*x$dp[x$g==PANEL[nm]]/med, 2)
  d}))
g  <- fread(file.path(C,"erg11_genotypes.csv"))
ph <- fread(file.path(B,"doc/Candida/tables/ctropicalis_fluconazole_phenotype.csv"))
m <- cn %>% inner_join(g, by="accession") %>%
     inner_join(ph[,.(accession=strain_name, mic_used, log2_mic, R_call, censored,
                      published_genotype=resistance_genotype, country, biosample, sra_run)], by="accession") %>%
     filter(status=="ok") %>%
     mutate(qc_pass = median_depth >= 20,
            ERG11_amplified = CN_ERG11 > 2.5,
            mechanism = case_when(Y132F==1 & ERG11_amplified ~ "Y132F + amplification",
                                  Y132F==1 ~ "Y132F only",
                                  ERG11_amplified ~ "amplification only",
                                  TRUE ~ "neither"))
fwrite(m, file.path(out,"tables/ctropicalis_analysis_table.csv"))
a <- m %>% filter(qc_pass)
cat("=== ANALYSIS SET\n"); cat("processed:", nrow(m), " passing QC (>=20x):", nrow(a), "\n")
print(table(a$R_call))
cat("\n=== mechanism vs resistance\n")
print(a %>% group_by(mechanism) %>% summarise(n=n(), R=sum(R_call=="R"), S=sum(R_call=="S"),
      median_MIC=median(mic_used), .groups="drop") %>% arrange(desc(median_MIC)) %>% as.data.frame())
cat("\nvariance in log2 MIC explained by the two ERG11 mechanisms:",
    round(100*summary(lm(log2_mic ~ Y132F + ERG11_amplified, data=a))$r.squared,1), "%\n")
cat("Y132F alone:", round(100*summary(lm(log2_mic ~ Y132F, data=a))$r.squared,1),
    "%   amplification alone:", round(100*summary(lm(log2_mic ~ ERG11_amplified, data=a))$r.squared,1), "%\n")
cat("\nresistant isolates with NEITHER mechanism:", sum(a$R_call=="R" & a$mechanism=="neither"),
    "of", sum(a$R_call=="R"), "\n")
