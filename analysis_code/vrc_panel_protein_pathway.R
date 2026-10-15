#!/usr/bin/env Rscript
# Link the literature-verified azole-resistance gene regions to protein families and pathways.
# Protein/pathway assignments come from eggNOG-mapper run on the Af293 proteins (not from memory);
# literature evidence columns come from doc/azole_resistance_gene_panel.md.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; P <- file.path(base,"analysis/vrc_gwas/05_panel_annotation")
ann <- read.delim(file.path(P,"panel.emapper.annotations"), comment.char="#", header=FALSE, stringsAsFactors=FALSE)
hdr <- c("query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description",
         "Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass",
         "BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs")
names(ann) <- hdr[seq_len(ncol(ann))]
ann <- ann %>% separate(query, c("gene","locus_tag","protein_id"), sep="\\|", remove=FALSE)
loci <- read.delim(file.path(P,"panel_loci.tsv"), stringsAsFactors=FALSE)

# curated layer: mechanism tier + literature evidence (from the verified panel reference)
lit <- tribble(
 ~gene,~rank,~mechanism_tier,~protein_group,~vrc_effect,~evidence,~key_variants,~citation,
 "cyp51A",1,"Drug target","Cytochrome P450, sterol 14a-demethylase","TR46/Y121F/T289A: MIC >=16 mg/L in 100% of 98 isolates","(a) validated","TR34/TR46/TR53/TR120 promoter; L98H, Y121F, T289A, G448S, G54, G138, M220","Song 2025 Lancet Microbe PMID 40651492",
 "hmg1",2,"Sterol pathway (upstream)","HMG-CoA reductase, sterol-sensing domain","4-8x across triazoles; bidirectional CRISPR","(a) validated","F262del, S305P, I412S","Rybak 2019 mBio PMID 30940706",
 "cdr1B",3,"Efflux","ABC transporter, PDR/CDR subfamily","knockout drops VRC 1 -> 0.25 mg/L (4x)","(a) validated","R339H (2x); largely expression-level","Fraczek 2013 JAC PMID 23580559",
 "atrR",4,"Transcriptional regulator","Zn2Cys6 fungal TF","dAtrR drops VRC 0.5 -> <0.015 mg/L","(a) validated, no clinical variant","none published - prospecting","Hagiwara 2017 PLoS Pathog PMID 28052140",
 "hapE",5,"Transcriptional regulator","CCAAT-binding complex, histone-fold","2x (0.25 -> 1 mg/L)","(a) validated, rare","P88L (2 isolates ever)","Camps 2012 PLoS ONE 10.1371/journal.pone.0050034",
 "srbA",6,"Transcriptional regulator","SREBP, basic helix-loop-helix","dSrbA drops VRC 0.5 -> 0.03; induction raises 8x","(a) validated","L250P (background-dependent)","Gsaller 2016 PLoS Pathog PMID 27438727",
 "ptaB",7,"Transcriptional regulator","LIM-binding transcriptional cofactor","4x VRC/ITC/POS/ISA on loss of function","(a) validated, single study","Q312stop, Q397fs, deletion","Handelman 2026 AAC PMID 41778915",
 "cyp51B",8,"Drug target (paralog)","Cytochrome P450, sterol 14a-demethylase","G457S: 2-4x VRC on reintroduction","(a) validated, low yield","G457S","Handelman 2026 AAC PMID 41778915",
 "erg25B",9,"Sterol pathway","C-4 methyl sterol oxidase (fatty acid hydroxylase family)","2-4x VRC/ISA/ITC/POS","(a) validated, single study","D253G, G178D, K294fs","Handelman 2026 AAC PMID 41778915",
 "rttA",10,"Transcriptional regulator","Zn2Cys6 fungal TF (Upc2-like)","overexpression 8x VRC; deletion 2x lower","(a) validated, no clinical variant","none published; REFERENCE CDS TRUNCATED","Birstonas 2025 mBio PMID 40937856",
 "nctA",11,"Transcriptional regulator","Negative cofactor 2, histone-fold","azoles + amphotericin B + terbinafine (figure-only)","(a) validated, MIC values unverified","null mutants","Furukawa 2020 Nat Commun PMID 31969561",
 "nctB",11,"Transcriptional regulator","Negative cofactor 2, histone-fold","as nctA","(a) validated, MIC values unverified","null mutants","Furukawa 2020 Nat Commun PMID 31969561",
 "complexI_29.9kDa",12,"Respiration","NADH dehydrogenase (complex I) subunit","deletion: ITC 0.25 -> >8; 4x VRC","(b) association + knockout","E180D in 12/22 R vs 2/20 S (p=0.0028)","Bromley 2016 PLoS One PMID 27438017",
 "cox10",13,"Respiration","Haem A farnesyltransferase (UbiA family)","ITC-centred; reduces drug uptake","(a) lab-selected allele","R243Q (lab); P17S, A423V (clinical, association only)","see panel reference",
 "rbdB",14,"Membrane proteolysis","Rhomboid intramembrane serine protease","dRbdB VRC 0.008 vs 0.19-0.25 (hypersusceptible)","(a) validated, no clinical variant","none","see panel reference",
 "mshA",15,"DNA mismatch repair","MutS family","MICs UNCHANGED - hypermutator, not an MIC predictor","(a) validated as emergence locus","not applicable","dos Reis 2019; Bottery 2024",
 "mdr1",16,"Efflux","ABC transporter, MDR/TAP subfamily","no VRC effect (2x ITC/POS only)","(a) validated, negative for VRC","-","see panel reference")

kegg_names <- c(ko00100="Steroid biosynthesis", ko00900="Terpenoid backbone biosynthesis",
  ko02010="ABC transporters", ko00190="Oxidative phosphorylation", ko00860="Porphyrin metabolism",
  ko03430="Mismatch repair", ko04612="Antigen processing (mammalian; NF-Y homology)",
  ko01100="Metabolic pathways", ko01110="Biosynthesis of secondary metabolites", ko01130="Antibiotic biosynthesis")
tab <- lit %>% left_join(ann %>% select(gene, locus_tag, protein_id, COG_category, PFAMs, KEGG_ko, KEGG_Pathway, EC, Description), by="gene") %>%
  left_join(loci %>% select(gene, panel_note=note), by="gene") %>%
  mutate(kegg_core = sapply(strsplit(KEGG_Pathway, ","), function(x){
           x <- x[grepl("^ko", x)]; x <- x[x %in% names(kegg_names)]
           x <- x[!x %in% c("ko01100","ko01110","ko01130")]
           if(!length(x)) "-" else paste(paste0(kegg_names[x]," (",x,")"), collapse="; ")})) %>%
  arrange(rank, gene) %>%
  select(rank, gene, locus_tag, protein_id, mechanism_tier, protein_group, pfam_domains=PFAMs, COG=COG_category,
         kegg_ko=KEGG_ko, pathway=kegg_core, vrc_effect, evidence, key_variants, citation, panel_note)
write.csv(tab, file.path(base,"output/table/azole_resistance_panel_annotation.csv"), row.names=FALSE, na="")
cat("=== panel annotation, ", nrow(tab), " regions\n", sep="")
print(as.data.frame(tab %>% select(rank, gene, mechanism_tier, protein_group, pathway)), right=FALSE)
cat("\n=== protein families represented:\n"); print(sort(table(tab$protein_group), decreasing=TRUE))
cat("\n=== mechanism tiers:\n"); print(table(tab$mechanism_tier))

# ---- figure: gene -> protein family -> pathway tier map
pd <- tab %>% filter(rank <= 13) %>%
  mutate(tier = factor(mechanism_tier, levels=c("Drug target","Drug target (paralog)","Sterol pathway (upstream)",
         "Sterol pathway","Efflux","Transcriptional regulator","Respiration")),
         path_short = case_when(grepl("Steroid", pathway) ~ "Ergosterol / steroid biosynthesis (ko00100)",
                                grepl("Terpenoid", pathway) ~ "Terpenoid backbone / mevalonate (ko00900)",
                                grepl("ABC", pathway) ~ "ABC transporters (ko02010)",
                                grepl("Oxidative", pathway) ~ "Oxidative phosphorylation (ko00190)",
                                TRUE ~ "No KEGG map (fungal transcriptional regulation)"),
         lab = paste0(gene, "  ", ifelse(rank<=10, paste0("#",rank), "")))
p <- ggplot(pd, aes(x=tier, y=reorder(lab, -rank), fill=path_short)) +
  geom_tile(colour="white", linewidth=1.1, width=.92, height=.86) +
  geom_text(aes(label=gsub(",.*","",pfam_domains)), size=2.5, colour="grey15", vjust=2.4) +
  scale_fill_manual(values=c("Ergosterol / steroid biosynthesis (ko00100)"="#7FB3A4",
    "Terpenoid backbone / mevalonate (ko00900)"="#B8D3C6","ABC transporters (ko02010)"="#E8C48A",
    "Oxidative phosphorylation (ko00190)"="#C9A3C4","No KEGG map (fungal transcriptional regulation)"="#A8BFD4")) +
  labs(x=NULL, y=NULL, fill="KEGG pathway",
       title="Azole-resistance gene regions: mechanism tier, protein family and pathway",
       subtitle="13 literature-verified regions; tile label = leading Pfam domain (eggNOG-mapper on Af293 proteins)") +
  theme_bw(base_size=10) + theme(axis.text.x=element_text(angle=20, hjust=1), legend.position="bottom",
    legend.direction="vertical", panel.grid=element_blank())
ggsave(file.path(base,"output/figure/azole_resistance_panel_pathways.pdf"), p, width=10, height=7.5)
cat("\nwrote output/table/azole_resistance_panel_annotation.csv and output/figure/azole_resistance_panel_pathways.pdf\n")
