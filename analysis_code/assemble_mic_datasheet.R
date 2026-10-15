#!/usr/bin/env Rscript
# Merge agent-extracted MIC tables into one data sheet (PASS-filtered genomes only).
# Outputs: output/table/aspergillus_mic_long.csv, aspergillus_mic_wide.csv, aspergillus_mic_sources.csv,
#          aspergillus_mic_coverage.csv and output/table/Aspergillus_MIC_datasheet.xlsx
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(openxlsx)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
src  <- "/tmp/claude-1000/-home-huanglab/25fa2051-9b57-45bc-894f-755b7f717add/scratchpad/mic"
files <- c(barber2021="barber2021.csv", fumigatus_clinical="fumigatus_clinical.csv", asia_cryptic="asia_cryptic.csv",
           flavus_niger_mixed="flavus_niger_mixed.csv", literature_sweep="literature_sweep.csv")
rd <- function(f) read.csv(file.path(src,f), stringsAsFactors=FALSE, na.strings=c("","NA"), colClasses="character")
d <- bind_rows(lapply(names(files), function(k) rd(files[k]) %>% mutate(extract_file=k)))
# --- source keys: the barber file used sample names as source_id -> normalise; drop sweep's duplicate of Barber 2021
d <- d %>% mutate(source_key = case_when(
        extract_file=="barber2021" ~ "barber2021",
        grepl("^PRJ", source_id) ~ source_id,
        TRUE ~ source_id)) %>%
     filter(!(extract_file=="literature_sweep" & source_key=="barber2021"))
# same paper extracted by two agents under different keys -> one key; keep the dedicated extraction over the sweep copy
alias <- c(brazil2026="PRJNA1354845", steenwyk2021="PRJNA673120", wu2020eid="PRJDB20599", PRJDB20599="PRJDB20599", drott2022="PMID36318036", PMID36318036="PMID36318036", PMID33080784="PRJNA836245", PRJNA836245="PRJNA836245", PRJNA1354845="PRJNA1354845", PRJNA673120="PRJNA673120")
d <- d %>% mutate(source_key = ifelse(source_key %in% names(alias), alias[source_key], source_key),
                  prio = ifelse(extract_file=="literature_sweep", 2, 1)) %>%
     group_by(accession, drug, source_key) %>% filter(prio==min(prio)) %>% ungroup() %>% select(-prio)
# --- PASS genomes only, and accession must physically exist in raw_data/coverage_pass
pass_acc <- basename(list.dirs(file.path(base,"raw_data/coverage_pass"), recursive=TRUE))
pass_acc <- pass_acc[grepl("^GC[AF]_", pass_acc)]
d <- d %>% filter(status=="PASS", accession %in% pass_acc)
# GCA_012897115.1 ("NRRL 2999") is a misidentified assembly slated for withdrawal (PMID 34854709); the published MIC belongs to a different isolate
d <- d %>% filter(accession != "GCA_012897115.1")
# --- normalise fields
d <- d %>% mutate(drug = tolower(trimws(drug)),
                  drug = recode(drug, "amphotericin b"="amphotericin_B", "amphotericin_b"="amphotericin_B"),
                  drug = gsub(" ", "_", drug),
                  mic_mg_per_L = suppressWarnings(as.numeric(mic_mg_per_L)),
                  mic_operator = ifelse(is.na(mic_operator) & !is.na(mic_mg_per_L), "=", mic_operator),
                  mic_text = ifelse(is.na(mic_mg_per_L), NA, paste0(ifelse(mic_operator=="=","",mic_operator), mic_mg_per_L)),
                  log2_mic = ifelse(is.na(mic_mg_per_L), NA, log2(mic_mg_per_L)),
                  interpretation = toupper(interpretation)) %>%
     distinct(accession, drug, source_key, mic_operator, mic_mg_per_L, interpretation, .keep_all=TRUE) %>%
     filter(!(is.na(mic_mg_per_L) & is.na(interpretation)))
# --- attach genome metadata from filter table
meta <- read.csv(file.path(base,"output/table/assembly_coverage_filter.csv"), stringsAsFactors=FALSE) %>%
        select(accession, genome_species=species, genome_strain=strain, level, coverage_x, contig_n50, bioproject, isolation_source, host, geo)
long <- d %>% rename(species_in_source=species) %>% left_join(meta, by="accession") %>% mutate(species=genome_species) %>%
  select(species, species_in_source, accession, strain_name, genome_strain, biosample, bioproject, drug, mic_operator, mic_mg_per_L, mic_text, log2_mic,
         interpretation, method, cyp51a_genotype, isolation_type, country, year, source_key, source_citation, source_doi_url, supplement_url,
         notes, level, coverage_x, contig_n50, isolation_source, host, geo) %>%
  arrange(species, accession, drug, source_key)
write.csv(long, file.path(base,"output/table/aspergillus_mic_long.csv"), row.names=FALSE, na="")
# --- wide: one row per genome; per drug the MIC (collapse multiple sources with " | ") and R/S
wide_mic <- long %>% filter(!is.na(mic_text)) %>% group_by(accession, drug) %>%
  summarise(v=paste(unique(mic_text), collapse=" | "), .groups="drop") %>% pivot_wider(names_from=drug, values_from=v, names_prefix="MIC_")
wide_rs  <- long %>% filter(!is.na(interpretation)) %>% group_by(accession, drug) %>%
  summarise(v=paste(unique(interpretation), collapse=" | "), .groups="drop") %>% pivot_wider(names_from=drug, values_from=v, names_prefix="RS_")
strain_info <- long %>% group_by(accession) %>% summarise(
  species=first(species), species_in_source=paste(unique(na.omit(species_in_source)), collapse=" | "), strain_name=paste(unique(na.omit(strain_name)), collapse=" | "), biosample=first(biosample), bioproject=first(bioproject),
  cyp51a_genotype=paste(unique(na.omit(cyp51a_genotype)), collapse=" | "), isolation_type=paste(unique(na.omit(isolation_type)), collapse=" | "),
  country=paste(unique(na.omit(country)), collapse=" | "), year=paste(unique(na.omit(year)), collapse=" | "),
  method=paste(unique(na.omit(method)), collapse=" | "), n_sources=n_distinct(source_key), sources=paste(unique(source_key), collapse="; "),
  source_doi_url=paste(unique(na.omit(source_doi_url)), collapse=" ; "), level=first(level), coverage_x=first(coverage_x), .groups="drop")
drug_order <- c("itraconazole","voriconazole","posaconazole","isavuconazole","amphotericin_B","caspofungin","micafungin","anidulafungin",
                "terbinafine","tebuconazole","difenoconazole","prochloraz","epoxiconazole","metconazole","propiconazole","natamycin","fluconazole","5-flucytosine")
ord <- function(cols, pre) { dr <- sub(paste0("^",pre),"",cols); cols[order(match(dr, drug_order, nomatch=99), dr)] }
wide <- strain_info %>% left_join(wide_mic, by="accession") %>% left_join(wide_rs, by="accession")
wide <- wide[, c(names(strain_info), ord(grep("^MIC_",names(wide),value=TRUE),"MIC_"), ord(grep("^RS_",names(wide),value=TRUE),"RS_"))] %>% arrange(species, accession)
write.csv(wide, file.path(base,"output/table/aspergillus_mic_wide.csv"), row.names=FALSE, na="")
# --- sources
sources <- long %>% group_by(source_key) %>% summarise(
  species=paste(unique(species), collapse=";"), n_genomes=n_distinct(accession), n_rows=n(),
  n_numeric_MIC=sum(!is.na(mic_mg_per_L)), n_RS_calls=sum(!is.na(interpretation)),
  drugs=paste(sort(unique(drug)), collapse="; "), method=paste(unique(na.omit(method)), collapse="; "),
  citation=first(na.omit(source_citation)), doi_url=first(na.omit(source_doi_url)), supplement_url=paste(unique(na.omit(supplement_url)), collapse=" ; "), .groups="drop") %>%
  arrange(desc(n_genomes))
write.csv(sources, file.path(base,"output/table/aspergillus_mic_sources.csv"), row.names=FALSE, na="")
# --- coverage: PASS genomes per species with/without MIC
allpass <- read.csv(file.path(base,"output/table/assembly_coverage_filter.csv"), stringsAsFactors=FALSE) %>% filter(status=="PASS")
cov <- allpass %>% group_by(species) %>% summarise(PASS_genomes=n(), .groups="drop") %>%
  left_join(long %>% group_by(species) %>% summarise(with_any_MIC_or_RS=n_distinct(accession),
            with_numeric_MIC=n_distinct(accession[!is.na(mic_mg_per_L)]), RS_call_only=n_distinct(accession[is.na(mic_mg_per_L)]) - 0, .groups="drop"), by="species") %>%
  mutate(across(where(is.numeric), ~replace_na(.x,0)), without_MIC = PASS_genomes - with_any_MIC_or_RS)
# RS_call_only should count genomes with NO numeric MIC at all
rs_only <- long %>% group_by(species, accession) %>% summarise(num=any(!is.na(mic_mg_per_L)), .groups="drop") %>% group_by(species) %>% summarise(RS_call_only=sum(!num), .groups="drop")
cov <- cov %>% select(-RS_call_only) %>% left_join(rs_only, by="species") %>% mutate(RS_call_only=replace_na(RS_call_only,0))
write.csv(cov, file.path(base,"output/table/aspergillus_mic_coverage.csv"), row.names=FALSE)
# --- per drug × species counts
perdrug <- long %>% group_by(species, drug) %>% summarise(n_genomes=n_distinct(accession), n_numeric=n_distinct(accession[!is.na(mic_mg_per_L)]),
           n_R=n_distinct(accession[interpretation %in% "R"]), n_S=n_distinct(accession[interpretation %in% "S"]), .groups="drop") %>% arrange(species, desc(n_genomes))
# --- Excel
readme <- data.frame(item=c("Scope","Rows","MIC units","mic_operator","interpretation","Sheets","Caveats","Generated by"),
  detail=c("Antifungal MICs from published sources for Aspergillus genomes held in raw_data/coverage_pass/ (assembly-quality PASS only). One genome may have several sources.",
           "MIC_long: one row per genome x drug x source. MIC_wide: one row per genome, MIC_<drug> = MIC in mg/L (multiple sources separated by ' | '), RS_<drug> = R/S call given by the source paper.",
           "mg/L (= ug/mL). log2_mic provided. Echinocandin values are MECs where the source used MEC (see notes).",
           "= exact; > or >= above highest tested; <= at/below lowest tested.",
           "R/S/I exactly as stated by the source; blank when the paper reports numbers only. Barber 2021 gives R/S only (numeric MICs only for its resistant isolates via Barber 2020 mBio).",
           "README, MIC_wide, MIC_long, Per_drug_counts, Sources, Genome_coverage.",
           "Methods differ between sources (EUCAST vs CLSI vs other) - see method column; harmonise before pooling. Rows flagged in notes as paywalled/unverified should be checked against the publisher PDF.",
           paste("analysis_code/assemble_mic_datasheet.R,", format(Sys.time(), "%Y-%m-%d %H:%M"))))
wb <- createWorkbook()
add <- function(name, df) { addWorksheet(wb, name); writeData(wb, name, df, headerStyle=createStyle(textDecoration="bold", fgFill="#DDEBF7")); freezePane(wb, name, firstRow=TRUE); setColWidths(wb, name, cols=1:ncol(df), widths="auto") }
add("README", readme); add("MIC_wide", wide); add("MIC_long", long); add("Per_drug_counts", perdrug); add("Sources", sources); add("Genome_coverage", cov)
saveWorkbook(wb, file.path(base,"output/table/Aspergillus_MIC_datasheet.xlsx"), overwrite=TRUE)
cat("long rows:", nrow(long), " genomes:", n_distinct(long$accession), " sources:", nrow(sources), "\n"); print(as.data.frame(cov)); print(as.data.frame(sources[,c("source_key","species","n_genomes","n_numeric_MIC","n_RS_calls")]))
