#!/usr/bin/env Rscript
# Phenotype tables from NUMERIC MIC values only (categorical R/S-only records are excluded).
# Headline drug: voriconazole; itraconazole/posaconazole/amphotericin B built the same way.
# Censoring rule (user-specified): a reported ">x" is analysed as 2x (one doubling dilution above the
#   highest concentration tested). The mirrored convention is applied to "<=x" -> x/2. The reported
#   value and operator are kept in mic_reported / mic_operator; mic_used is the analysis value.
# Continuous phenotype  : log2(mic_used in mg/L).
# Binary phenotype      : EUCAST v10 A. fumigatus clinical breakpoints, R if mic_used > breakpoint
#                         (VRC >1, ITC >1, POS >0.25, AMB >1). Alternative >2 call kept as R_strict.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(openxlsx)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
long <- read.csv(file.path(base,"output/table/aspergillus_mic_long.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
num  <- long %>% filter(!is.na(mic_mg_per_L)) %>%                # <<< numeric only
  mutate(mic_reported = mic_mg_per_L,
         mic_used = case_when(mic_operator %in% c(">",">=") ~ mic_mg_per_L*2,
                              mic_operator == "<="          ~ mic_mg_per_L/2,
                              TRUE                          ~ mic_mg_per_L),
         censor_rule = case_when(mic_operator %in% c(">",">=") ~ "censored high: reported value doubled",
                                 mic_operator == "<="          ~ "censored low: reported value halved",
                                 TRUE                          ~ ""))
write.csv(num, file.path(base,"output/table/aspergillus_mic_numeric_only.csv"), row.names=FALSE, na="")
secondary <- "celia2024"                                          # compiled meta-table -> lowest priority
bp <- c(voriconazole=1, itraconazole=1, posaconazole=0.25, amphotericin_B=1)

collapse_one <- function(df) {                                    # one row per genome x drug
  df %>% mutate(prio=ifelse(source_key %in% secondary, 2, 1)) %>%
    arrange(accession, prio, desc(mic_used)) %>% group_by(accession) %>%
    summarise(mic_used=first(mic_used), mic_reported=first(mic_reported), mic_operator=first(mic_operator),
              censor_rule=first(censor_rule), mic_source=first(source_key),
              mic_method=first(method), n_sources=n_distinct(source_key),
              mic_all_sources=paste(unique(paste0(source_key,":",ifelse(mic_operator=="=","",mic_operator),mic_reported)), collapse=" | "), .groups="drop")
}
meta <- num %>% group_by(accession) %>% summarise(
  species=first(species), strain_name=paste(unique(na.omit(strain_name)), collapse=" | "), biosample=first(biosample), bioproject=first(bioproject),
  cyp51a_genotype=paste(unique(na.omit(cyp51a_genotype)), collapse=" | "), isolation_type=paste(unique(na.omit(isolation_type)), collapse=" | "),
  country=paste(unique(na.omit(country)), collapse=" | "), year=paste(unique(na.omit(year)), collapse=" | "),
  level=first(level), coverage_x=first(coverage_x), contig_n50=first(contig_n50),
  source_doi_url=paste(unique(na.omit(source_doi_url)), collapse=" ; "), .groups="drop")

build <- function(dr, sp="fumigatus") {
  d <- num %>% filter(drug==dr, species==sp)
  if (nrow(d)==0) return(NULL)
  b <- bp[[dr]]
  collapse_one(d) %>% left_join(meta, by="accession") %>%
    mutate(drug=dr, breakpoint_mg_per_L=b, censored=mic_operator %in% c(">",">=","<="),
           log2_mic=log2(mic_used),
           R_call=ifelse(mic_used>b, "R","S"),
           R_strict=ifelse(mic_used>2, "R","S"),
           pheno_binary=as.integer(R_call=="R"),
           genome_dir=file.path(base,"raw_data/coverage_pass",paste0("Aspergillus_",sp),accession)) %>%
    select(accession, species, strain_name, biosample, bioproject, drug, mic_used, mic_reported, mic_operator, censored, censor_rule, log2_mic,
           breakpoint_mg_per_L, R_call, pheno_binary, R_strict, mic_source, mic_method, n_sources, mic_all_sources,
           cyp51a_genotype, isolation_type, country, year, level, coverage_x, contig_n50, source_doi_url, genome_dir) %>%
    arrange(desc(pheno_binary), mic_used)
}
drugs <- c("voriconazole","itraconazole","posaconazole","amphotericin_B")
tabs <- setNames(lapply(drugs, build), drugs)
vrc <- tabs$voriconazole
write.csv(vrc, file.path(base,"output/table/voriconazole_phenotype_fumigatus.csv"), row.names=FALSE, na="")
write.table(vrc %>% select(samples=accession, voriconazole_R=pheno_binary), file.path(base,"output/table/voriconazole_pheno_binary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(vrc %>% select(samples=accession, voriconazole_log2MIC=log2_mic), file.path(base,"output/table/voriconazole_pheno_log2mic.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
allaz <- bind_rows(tabs)
write.csv(allaz, file.path(base,"output/table/azole_phenotypes_numeric_fumigatus.csv"), row.names=FALSE, na="")
# combined multi-drug matrix (one row per genome)
wide <- allaz %>% select(accession, drug, mic_used) %>% pivot_wider(names_from=drug, values_from=mic_used, names_prefix="MIC_") %>%
  left_join(allaz %>% select(accession, drug, R_call) %>% pivot_wider(names_from=drug, values_from=R_call, names_prefix="R_"), by="accession") %>%
  left_join(meta %>% select(accession, strain_name, cyp51a_genotype, country, year, isolation_type), by="accession")
write.csv(wide, file.path(base,"output/table/azole_phenotypes_numeric_wide.csv"), row.names=FALSE, na="")
# ---- report
cat("NUMERIC-ONLY phenotypes, A. fumigatus\n")
for (dr in drugs) { t <- tabs[[dr]]; cat(sprintf("%-15s n=%3d  R=%2d  S=%3d  (breakpoint >%s mg/L; censored %d)\n", dr, nrow(t), sum(t$pheno_binary), sum(!t$pheno_binary), bp[[dr]], sum(t$censored))) }
cat("\nVRC MIC distribution (analysis value, censored already doubled):\n"); print(table(vrc$mic_used))
cat("VRC censored values (reported -> used):\n"); print(vrc %>% filter(censored) %>% transmute(accession, strain_name, reported=paste0(mic_operator,mic_reported), used=mic_used) %>% as.data.frame())
cat("\nVRC: R vs cyp51A TR genotype:\n"); print(table(TR=grepl("TR3|TR4",vrc$cyp51a_genotype), R=vrc$R_call))
cat("\nVRC by country:\n"); print(table(ifelse(vrc$country=="","(unstated)",vrc$country)))
cat("\nVRC by method:\n"); print(table(ifelse(is.na(vrc$mic_method),"(unstated)",vrc$mic_method)))
cat("\nalternative breakpoint >2 mg/L: R =", sum(vrc$R_strict=="R"), "\n")
cat("\nother species with numeric MIC (not analysable, n too small):\n")
print(num %>% filter(species!="fumigatus") %>% group_by(species) %>% summarise(genomes=n_distinct(accession), drugs=paste(sort(unique(drug)),collapse=","), .groups="drop") %>% as.data.frame())
# ---- figure
pd <- vrc %>% mutate(lbl=factor(as.character(mic_used), levels=c("0.25","0.5","1","2","4","8","16","32","64")))
p1 <- ggplot(pd, aes(lbl, fill=R_call)) + geom_bar() + geom_vline(xintercept=3.5, linetype=2, colour="grey40") +
  scale_fill_manual(values=c(R="#b2182b", S="#2166ac")) +
  labs(x="Voriconazole MIC (mg/L)", y="Genomes", fill=NULL, title="Voriconazole, numeric MIC only",
       subtitle=paste0("n = ", nrow(vrc), " genomes; censored '>x' analysed as 2x; dashed line = EUCAST breakpoint (R > 1 mg/L)")) + theme_bw(base_size=11)
p2 <- ggplot(pd, aes(x=grepl("TR3|TR4",cyp51a_genotype), y=log2_mic, colour=R_call)) +
  geom_jitter(width=.15, height=.05, size=2, alpha=.8) + scale_colour_manual(values=c(R="#b2182b", S="#2166ac")) +
  scale_x_discrete(labels=c("FALSE"="no TR mutation","TRUE"="cyp51A TR34/TR46")) +
  labs(x=NULL, y="log2 MIC (mg/L)", colour=NULL, title="Positive control", subtitle="known resistance genotype vs measured MIC") + theme_bw(base_size=11)
pdf(file.path(base,"output/figure/voriconazole_phenotype_fumigatus.pdf"), width=11, height=4.5); gridExtra::grid.arrange(p1,p2,ncol=2); dev.off()
# ---- workbook sheets
xl <- file.path(base,"output/table/Aspergillus_MIC_datasheet.xlsx"); wb <- loadWorkbook(xl)
for (s in c("VRC_phenotype","Numeric_only_azoles")) if (s %in% names(wb)) removeWorksheet(wb, s)
addWorksheet(wb,"VRC_phenotype"); writeData(wb,"VRC_phenotype",vrc, headerStyle=createStyle(textDecoration="bold", fgFill="#FCE4D6")); freezePane(wb,"VRC_phenotype",firstRow=TRUE); setColWidths(wb,"VRC_phenotype",1:ncol(vrc),"auto")
addWorksheet(wb,"Numeric_only_azoles"); writeData(wb,"Numeric_only_azoles",wide, headerStyle=createStyle(textDecoration="bold", fgFill="#E2EFDA")); freezePane(wb,"Numeric_only_azoles",firstRow=TRUE); setColWidths(wb,"Numeric_only_azoles",1:ncol(wide),"auto")
saveWorkbook(wb, xl, overwrite=TRUE)
cat("\nwritten: voriconazole_phenotype_fumigatus.csv, azole_phenotypes_numeric_{fumigatus,wide}.csv, aspergillus_mic_numeric_only.csv, 2 tsv, figure, 2 xlsx sheets\n")
