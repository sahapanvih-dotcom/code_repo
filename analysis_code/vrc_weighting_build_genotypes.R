#!/usr/bin/env Rscript
# Build the gene-level genotype matrix for the weighting / variance-partition models.
# Three encodings per gene:
#   known  : literature-defined resistance alleles only (highest specificity, lowest sensitivity)
#   rare   : burden of non-synonymous substitutions at allele frequency <= 10% (filters the common
#            neutral haplotypes; stands in for domain-restricted burden, which needs Pfam HMMs that
#            are not available on this host)
#   burden : total non-synonymous substitution count (most sensitive, noisiest)
# Protein-length differences reported by the caller are treated as indel variants, not failures.
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; W <- file.path(base,"analysis/vrc_gwas")
out  <- file.path(base,"doc/vrc_gene_weighting_2026-08")

lt2gene <- c(AFUA_2G02690="atrR", AFUA_2G12910="ptaB", AFUA_7G04740="rttA", AFUA_2G14250="nctA",
             AFUA_3G02340="nctB", AFUA_2G10600="complexI", AFUA_4G08340="cox10")
p1 <- read.csv(file.path(W,"04_association/panel/panel_variants.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
p2 <- read.csv(file.path(W,"04_association/panel/panel_variants_extended.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")) %>%
      mutate(gene = unname(lt2gene[gene]))
pan <- bind_rows(p1, p2) %>% mutate(gene = recode(gene, abcC="cdr1B"))

# the 10 literature regions that carry variation here (+ the 4 invariant ones kept for reporting)
KEEP <- c("cyp51A","hmg1","cdr1B","srbA","cyp51B","atrR","ptaB","nctA","complexI","cox10")
INVAR <- c("hapE","erg25B","nctB","rttA")
pan <- pan %>% filter(gene %in% c(KEEP, INVAR))

# ---- indel flag from length differences
pan <- pan %>% mutate(indel = grepl("^length_diff", status),
                      callable = status == "ok" | indel)
cat("callable calls:", sum(pan$callable), "of", nrow(pan), "\n")

# ---- long form of substitutions
subs <- pan %>% filter(status=="ok", !is.na(substitutions)) %>%
  separate_rows(substitutions, sep=";") %>% filter(substitutions != "") %>%
  transmute(accession, gene, sub = substitutions)
freq <- subs %>% count(gene, sub, name="carriers") %>%
  mutate(af = carriers / n_distinct(pan$accession))
write.csv(freq %>% arrange(gene, desc(carriers)), file.path(out,"tables/substitution_frequencies.csv"), row.names=FALSE)
cat("distinct substitutions:", nrow(freq), " rare (AF<=0.10):", sum(freq$af<=0.10), "\n")

# ---- encoding: literature-known resistance alleles
cy <- read.csv(file.path(base,"output/table/vrc_cyp51A_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
KNOWN <- list(
  cyp51A = c("L98H","Y121F","T289A","G448S","G54E","G54R","G54V","G54W","G138C","G138S","M220I","M220K","M220T","M220V","P216L","F219I","F219C"),
  hmg1   = c("S305P","I412S","F262del"),
  cdr1B  = c("R339H"),
  cyp51B = c("G457S"),
  srbA   = c("L250P"),
  complexI = c("E180D"),
  erg25B = c("D253G","G178D"))
known_tab <- subs %>% mutate(is_known = mapply(function(g,s) s %in% (if (is.null(KNOWN[[g]])) character(0) else KNOWN[[g]]), gene, sub))

accs <- sort(unique(pan$accession))
mk <- function(df, valcol) df %>% complete(accession=accs, gene=KEEP, fill=setNames(list(0), valcol))

enc_known <- known_tab %>% filter(is_known) %>% count(accession, gene, name="v") %>%
  bind_rows(tibble(accession=accs[1], gene=KEEP, v=0)) %>% group_by(accession, gene) %>%
  summarise(v=max(v), .groups="drop") %>% complete(accession=accs, gene=KEEP, fill=list(v=0)) %>%
  mutate(v = as.integer(v>0))
# cyp51A promoter TR is a known mechanism but is not a substitution -> add it in
tr <- cy %>% transmute(accession, tr = as.integer(promoter_genotype != "none"))
enc_known <- enc_known %>% left_join(tr, by="accession") %>%
  mutate(v = ifelse(gene=="cyp51A", pmax(v, ifelse(is.na(tr),0,tr)), v)) %>% select(-tr)

enc_rare <- subs %>% left_join(freq %>% select(gene, sub, af), by=c("gene","sub")) %>%
  filter(af <= 0.10) %>% count(accession, gene, name="v") %>%
  complete(accession=accs, gene=KEEP, fill=list(v=0)) %>% filter(gene %in% KEEP)
enc_burden <- subs %>% count(accession, gene, name="v") %>%
  complete(accession=accs, gene=KEEP, fill=list(v=0)) %>% filter(gene %in% KEEP)
# indels count toward rare and burden
ind <- pan %>% filter(indel, gene %in% KEEP) %>% count(accession, gene, name="iv")
add_ind <- function(df) df %>% left_join(ind, by=c("accession","gene")) %>%
  mutate(v = v + ifelse(is.na(iv),0,iv)) %>% select(-iv)
enc_rare <- add_ind(enc_rare); enc_burden <- add_ind(enc_burden)

wide <- function(df, pre) df %>% filter(gene %in% KEEP) %>%
  pivot_wider(names_from=gene, values_from=v, names_prefix=pre, values_fill=0) %>% arrange(accession)
G <- wide(enc_known,"known_") %>%
     left_join(wide(enc_rare,"rare_"), by="accession") %>%
     left_join(wide(enc_burden,"burden_"), by="accession")

ph <- read.csv(file.path(base,"output/table/voriconazole_phenotype_fumigatus.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")) %>%
      select(accession, strain_name, mic_used, log2_mic, R_call, mic_method, country)
qc <- read.csv(file.path(base,"output/table/vrc_phase0_qc.csv"), stringsAsFactors=FALSE) %>% select(accession, clone_group)
dat <- ph %>% left_join(qc, by="accession") %>% left_join(G, by="accession") %>%
       left_join(cy %>% select(accession, promoter_genotype, known_mechanism), by="accession")
write.csv(dat, file.path(out,"tables/model_input_genotypes.csv"), row.names=FALSE, na="")
cat("\nmodel input:", nrow(dat), "genomes x", ncol(dat), "columns\n")
cat("\ncarrier counts per gene (known / rare>0 / burden>0):\n")
for (g in KEEP) cat(sprintf("  %-10s %2d  %2d  %2d\n", g, sum(dat[[paste0("known_",g)]]>0),
    sum(dat[[paste0("rare_",g)]]>0), sum(dat[[paste0("burden_",g)]]>0)))
cat("\ninvariant regions (no testable variation):", paste(INVAR, collapse=", "), "\n")
