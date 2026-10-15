#!/usr/bin/env Rscript
# Build the voriconazole phenotype table for A. fumigatus PASS genomes (GWAS-ready).
# Primary design: binary R/S over all genomes with a VRC phenotype.
# Secondary: continuous log2(MIC) over the subset with a numeric MIC.
#
# Breakpoint policy (EUCAST v10.0+, A. fumigatus, voriconazole): S <= 1, R > 1, ATU = 2.
#   R_strict   : MIC > 2            (conservative; ATU counted as susceptible)
#   R_eucast   : MIC > 1            (current clinical breakpoint; ATU=2 counted as resistant)
#   R_final    : per-drug R/S as published where the source gives one (Barber 2021, EUCAST),
#                otherwise derived with R_eucast. This is the phenotype to analyse.
# NOTE: celia2024 is a compiled meta-table whose "resistance" label is isolate-level (identical
#   across itraconazole/voriconazole/posaconazole), NOT a per-drug call -> its interpretations are
#   discarded here; only its numeric MICs are used, and only when no primary source reports one.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(openxlsx)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
long <- read.csv(file.path(base,"output/table/aspergillus_mic_long.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
v <- long %>% filter(drug=="voriconazole", species=="fumigatus")

secondary <- c("celia2024")                      # compiled meta-tables (numeric MIC only, lowest priority)
num <- v %>% filter(!is.na(mic_mg_per_L)) %>%
  mutate(prio = ifelse(source_key %in% secondary, 2, 1)) %>%
  arrange(accession, prio, desc(mic_mg_per_L)) %>% group_by(accession) %>%
  summarise(mic_mg_per_L=first(mic_mg_per_L), mic_operator=first(mic_operator), mic_source=first(source_key),
            mic_method=first(method), mic_all=paste(unique(paste0(source_key,":",ifelse(mic_operator=="=","",mic_operator),mic_mg_per_L)), collapse=" | "),
            n_mic_sources=n_distinct(source_key), .groups="drop")
# per-drug published R/S (exclude the isolate-level flag from compiled tables)
rs <- v %>% filter(!is.na(interpretation), !source_key %in% secondary) %>%
  group_by(accession) %>% summarise(rs_published=paste(sort(unique(interpretation)), collapse="/"),
            rs_source=paste(unique(source_key), collapse=";"), rs_method=first(method), .groups="drop")
meta <- long %>% filter(species=="fumigatus") %>% group_by(accession) %>%
  summarise(strain_name=paste(unique(na.omit(strain_name)), collapse=" | "), biosample=first(biosample), bioproject=first(bioproject),
            cyp51a_genotype=paste(unique(na.omit(cyp51a_genotype)), collapse=" | "),
            isolation_type=paste(unique(na.omit(isolation_type)), collapse=" | "), country=paste(unique(na.omit(country)), collapse=" | "),
            year=paste(unique(na.omit(year)), collapse=" | "), level=first(level), coverage_x=first(coverage_x), contig_n50=first(contig_n50),
            sources=paste(sort(unique(source_key)), collapse=";"), source_doi_url=paste(unique(na.omit(source_doi_url)), collapse=" ; "), .groups="drop")

ph <- v %>% distinct(accession) %>% left_join(num, by="accession") %>% left_join(rs, by="accession") %>% left_join(meta, by="accession") %>%
  mutate(censored = mic_operator %in% c(">",">=","<="),
         log2_mic = ifelse(is.na(mic_mg_per_L), NA, log2(mic_mg_per_L)),
         R_strict = case_when(is.na(mic_mg_per_L) ~ NA_character_,
                              mic_operator %in% c(">",">=") & mic_mg_per_L>=2 ~ "R",
                              mic_mg_per_L > 2 ~ "R", TRUE ~ "S"),
         R_eucast = case_when(is.na(mic_mg_per_L) ~ NA_character_,
                              mic_operator %in% c(">",">=") & mic_mg_per_L>=1 ~ "R",
                              mic_mg_per_L > 1 ~ "R", TRUE ~ "S"),
         rs_published = ifelse(rs_published %in% c("R/S"), NA, rs_published),   # internally inconsistent -> treat as missing
         R_final = coalesce(rs_published, R_eucast),
         R_final_origin = ifelse(!is.na(rs_published), paste0("published (",rs_source,")"), ifelse(!is.na(R_eucast), "derived EUCAST >1 mg/L", NA)),
         pheno_binary = ifelse(is.na(R_final), NA_integer_, as.integer(R_final=="R")),
         genome_fasta = file.path(base,"raw_data/coverage_pass/Aspergillus_fumigatus",accession),
         agreement = case_when(is.na(rs_published)|is.na(R_eucast) ~ NA_character_, rs_published==R_eucast ~ "agree", TRUE ~ "DISAGREE")) %>%
  select(accession, strain_name, biosample, bioproject, pheno_binary, R_final, R_final_origin, rs_published, rs_source, rs_method,
         mic_mg_per_L, mic_operator, log2_mic, censored, R_eucast, R_strict, agreement, mic_source, mic_method, n_mic_sources, mic_all,
         cyp51a_genotype, isolation_type, country, year, level, coverage_x, contig_n50, sources, source_doi_url, genome_fasta) %>%
  arrange(desc(pheno_binary), accession)
write.csv(ph, file.path(base,"output/table/voriconazole_phenotype_fumigatus.csv"), row.names=FALSE, na="")
# GWAS-ready minimal files
write.table(ph %>% filter(!is.na(pheno_binary)) %>% select(samples=accession, voriconazole_R=pheno_binary),
            file.path(base,"output/table/voriconazole_pheno_binary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(ph %>% filter(!is.na(log2_mic)) %>% select(samples=accession, voriconazole_log2MIC=log2_mic),
            file.path(base,"output/table/voriconazole_pheno_log2mic.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
# ---- summaries
cat("genomes with VRC phenotype:", sum(!is.na(ph$pheno_binary)), " R:", sum(ph$pheno_binary==1,na.rm=TRUE), " S:", sum(ph$pheno_binary==0,na.rm=TRUE), "\n")
cat("with numeric MIC:", sum(!is.na(ph$mic_mg_per_L)), " censored:", sum(ph$censored, na.rm=TRUE), "\n")
cat("phenotype origin:\n"); print(table(ph$R_final_origin, useNA="ifany"))
cat("published vs derived agreement (where both exist):\n"); print(table(ph$agreement, useNA="ifany"))
cat("R under alternative breakpoints (numeric subset):  >1 mg/L:", sum(ph$R_eucast=="R",na.rm=TRUE), "  >2 mg/L:", sum(ph$R_strict=="R",na.rm=TRUE), "\n")
cat("duplicate strain names (possible same isolate):\n"); print(ph %>% filter(strain_name!="") %>% count(strain_name) %>% filter(n>1) %>% as.data.frame())
cat("cyp51A genotypes among R:\n"); print(table(ph$cyp51a_genotype[ph$pheno_binary==1 & ph$cyp51a_genotype!=""]))
cat("country:\n"); print(table(ifelse(ph$country=="","(unstated)",ph$country)))
# ---- figures
pd <- ph %>% filter(!is.na(mic_mg_per_L)) %>% mutate(lbl=factor(paste0(ifelse(mic_operator=="=","",mic_operator), mic_mg_per_L),
        levels=c("0.25","0.5","1","2","4","8",">8",">16","32",">32")))
p1 <- ggplot(pd, aes(lbl, fill=R_final)) + geom_bar() +
  geom_vline(xintercept=3.5, linetype=2, colour="grey40") +
  scale_fill_manual(values=c(R="#b2182b", S="#2166ac"), na.value="grey70") +
  labs(x="Voriconazole MIC (mg/L)", y="Genomes", fill="Phenotype",
       subtitle="dashed line = EUCAST clinical breakpoint (R > 1 mg/L)") + theme_bw(base_size=11)
p2 <- ph %>% filter(!is.na(pheno_binary)) %>% mutate(src=ifelse(grepl("^published", R_final_origin), "published R/S", "derived from MIC")) %>%
  ggplot(aes(src, fill=R_final)) + geom_bar() + scale_fill_manual(values=c(R="#b2182b", S="#2166ac")) +
  labs(x=NULL, y="Genomes", fill="Phenotype", subtitle=paste0("n = ", sum(!is.na(ph$pheno_binary)), " A. fumigatus genomes")) + theme_bw(base_size=11)
pdf(file.path(base,"output/figure/voriconazole_phenotype_fumigatus.pdf"), width=10, height=4.5); gridExtra::grid.arrange(p1,p2,ncol=2); dev.off()
# ---- append sheet to the datasheet workbook
xl <- file.path(base,"output/table/Aspergillus_MIC_datasheet.xlsx")
wb <- loadWorkbook(xl); if ("VRC_phenotype" %in% names(wb)) removeWorksheet(wb, "VRC_phenotype")
addWorksheet(wb, "VRC_phenotype"); writeData(wb, "VRC_phenotype", ph, headerStyle=createStyle(textDecoration="bold", fgFill="#FCE4D6"))
freezePane(wb, "VRC_phenotype", firstRow=TRUE); setColWidths(wb, "VRC_phenotype", cols=1:ncol(ph), widths="auto"); saveWorkbook(wb, xl, overwrite=TRUE)
cat("wrote voriconazole_phenotype_fumigatus.csv, pheno tsv x2, figure, and VRC_phenotype sheet\n")
