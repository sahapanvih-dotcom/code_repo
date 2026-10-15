#!/usr/bin/env Rscript
# Case series: full panel profile for every resistant genome that carries no known cyp51A mechanism.
# Descriptive, not inferential - no statistical power required, which is why it works at n=55.
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; W <- file.path(base,"analysis/vrc_gwas")
out <- file.path(base,"doc/vrc_gene_weighting_2026-08")
GENES <- c("cyp51A","cyp51B","hmg1","cdr1B","atrR","ptaB","nctA","srbA","complexI","cox10")
lt2 <- c(AFUA_2G02690="atrR", AFUA_2G12910="ptaB", AFUA_7G04740="rttA", AFUA_2G14250="nctA",
         AFUA_3G02340="nctB", AFUA_2G10600="complexI", AFUA_4G08340="cox10")
pan <- bind_rows(
  read.csv(file.path(W,"04_association/panel/panel_variants.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")),
  read.csv(file.path(W,"04_association/panel/panel_variants_extended.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")) %>%
    mutate(gene = unname(lt2[gene]))) %>% mutate(gene = recode(gene, abcC="cdr1B"))
cy <- read.csv(file.path(base,"output/table/vrc_cyp51A_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
ph <- read.csv(file.path(base,"output/table/voriconazole_phenotype_fumigatus.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
az <- read.csv(file.path(base,"output/table/azole_phenotypes_numeric_fumigatus.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
qc <- read.csv(file.path(base,"output/table/vrc_phase0_qc.csv"), stringsAsFactors=FALSE)
# allele frequencies, to mark which variants are rare
subs <- pan %>% filter(status=="ok", !is.na(substitutions)) %>% separate_rows(substitutions, sep=";") %>%
  filter(substitutions!="") %>% transmute(accession, gene, sub=substitutions)
freq <- subs %>% count(gene, sub, name="carriers") %>% mutate(af = carriers/55)

cases <- cy %>% filter(R_call=="R", !(known_mechanism %in% c("TRUE", TRUE))) %>% pull(accession)
cat("resistant genomes with no known cyp51A mechanism:", length(cases), "\n\n")
rows <- lapply(cases, function(a) {
  s <- subs %>% filter(accession==a) %>% left_join(freq, by=c("gene","sub")) %>%
       filter(gene %in% GENES, af <= 0.10) %>% arrange(gene, sub)
  ind <- pan %>% filter(accession==a, grepl("^length_diff", status), gene %in% GENES)
  info <- ph %>% filter(accession==a)
  oz <- az %>% filter(accession==a); others <- setNames(oz$mic_used, oz$drug)
  data.frame(accession=a, strain=info$strain_name, country=info$country, vrc_mic=info$mic_used,
    itc_mic=ifelse("itraconazole" %in% names(others), others[["itraconazole"]], NA),
    pos_mic=ifelse("posaconazole" %in% names(others), others[["posaconazole"]], NA),
    amb_mic=ifelse("amphotericin_B" %in% names(others), others[["amphotericin_B"]], NA),
    n_rare_variants=nrow(s),
    rare_variants=ifelse(nrow(s)==0, "none", paste(paste0(s$gene,":",s$sub," (n=",s$carriers,")"), collapse="; ")),
    indels=ifelse(nrow(ind)==0, "none", paste(paste0(ind$gene,":",ind$status), collapse="; ")),
    clone_group=qc$clone_group[match(a, qc$accession)], stringsAsFactors=FALSE)
}) %>% bind_rows() %>% arrange(desc(vrc_mic))
write.csv(rows, file.path(out,"tables/case_series_unexplained_resistant.csv"), row.names=FALSE, na="")
for (i in seq_len(nrow(rows))) { r <- rows[i,]
  cat(sprintf("--- %s  (%s, %s)\n    VRC %s | ITC %s | POS %s | AMB %s mg/L\n    rare panel variants (%d): %s\n    indels: %s\n\n",
      r$strain, r$accession, r$country, r$vrc_mic, r$itc_mic, r$pos_mic, r$amb_mic, r$n_rare_variants, r$rare_variants, r$indels))
}
cat("=== how many carry ANY rare panel variant:", sum(rows$n_rare_variants>0), "of", nrow(rows), "\n")
cat("=== genes recurring across cases:\n")
gg <- rows %>% filter(n_rare_variants>0) %>% separate_rows(rare_variants, sep="; ") %>%
  mutate(g=sub(":.*","",rare_variants)) %>% count(g, sort=TRUE)
print(as.data.frame(gg))
cat("\n=== cross-azole pattern (is resistance azole-wide or VRC-specific?):\n")
print(rows %>% transmute(strain, vrc_mic, itc_mic, pos_mic,
      pattern=ifelse(!is.na(itc_mic) & itc_mic>1 & !is.na(pos_mic) & pos_mic>0.25, "pan-azole",
              ifelse(!is.na(itc_mic) & itc_mic>1, "VRC+ITC", "VRC-predominant"))) %>% as.data.frame())
