suppressPackageStartupMessages({library(data.table); library(dplyr)})
C <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta/analysis/candida"
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
have <- sub("\\.vcf\\.gz$","",basename(Sys.glob(file.path(C,"vcf/*.vcf.gz"))))
fs <- file.path(C,"depth", paste0(have,".regions.bed.gz")); fs <- fs[file.exists(fs)]
d <- rbindlist(lapply(fs, function(f){
  x <- fread(f, col.names=c("c","s","e","g","dp")); med <- median(x$dp[x$dp>0])
  data.table(acc=sub("\\.regions\\.bed\\.gz$","",basename(f)), med=med,
             erg11=round(2*x$dp[x$g=="CTRG_05283"]/med,2))}))
d <- d[med>=20]
ph <- read.csv(file.path(B,"doc/Candida/tables/ctropicalis_fluconazole_phenotype.csv"), stringsAsFactors=FALSE)
ph$acc <- ph$strain_name
m <- merge(d, ph[,c("acc","mic_used","R_call","resistance_genotype")], by="acc")
pub <- m[grepl("copy number", m$resistance_genotype), ]
pub$pub_cn <- as.numeric(sub(".*copy number ([0-9]+).*", "\\1", pub$resistance_genotype))
pub <- pub[!is.na(pub$pub_cn), ]
cat("=== VALIDATION vs Fan 2023 published ERG11 copy number\n")
cat("isolates comparable:", nrow(pub), "\n")
cat("Pearson r =", round(cor(pub$erg11, pub$pub_cn),3),
    " Spearman rho =", round(cor(pub$erg11, pub$pub_cn, method="spearman"),3), "\n\n")
tb <- table(published_amplified = pub$pub_cn > 2, our_call_amplified = pub$erg11 > 2.5)
print(tb)
sens <- tb["TRUE","TRUE"]/sum(tb["TRUE",]); spec <- tb["FALSE","FALSE"]/sum(tb["FALSE",])
cat("\nsensitivity:", round(100*sens), "%   specificity:", round(100*spec), "%\n")
cat("\nour median CN by published CN:\n")
print(pub %>% group_by(pub_cn) %>% summarise(n=n(), our_median=round(median(erg11),2), .groups="drop") %>% as.data.frame())
