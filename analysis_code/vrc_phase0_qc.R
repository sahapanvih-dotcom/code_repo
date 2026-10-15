#!/usr/bin/env Rscript
# Phase 0 QC for the voriconazole GWAS set: species ID (fastANI vs references),
# clone/near-duplicate detection (all-vs-all ANI + mash), assembly stats.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; W <- file.path(base,"analysis/vrc_gwas")
nm <- function(x) sub("\\.fna$","",basename(x))
ph <- read.csv(file.path(base,"output/table/voriconazole_phenotype_fumigatus.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
# ---- species ID
ani_ref <- read.delim(file.path(W,"00_qc/ani_vs_refs.tsv"), header=FALSE,
            col.names=c("q","r","ani","frag_map","frag_tot")) %>% mutate(q=nm(q), r=sub("^ref_A_","",nm(r)))
spec <- ani_ref %>% group_by(q) %>% slice_max(ani, n=1) %>% ungroup() %>% select(accession=q, best_ref=r, ani_best=ani)
# ---- all-vs-all ANI
ava <- read.delim(file.path(W,"00_qc/ani_all_vs_all.tsv"), header=FALSE, col.names=c("q","r","ani","fm","ft")) %>%
       mutate(q=nm(q), r=nm(r)) %>% filter(q!=r)
# ---- mash
mash <- read.delim(file.path(W,"00_qc/mash_dist.tsv"), header=FALSE, col.names=c("q","r","dist","p","shared")) %>%
        mutate(q=nm(q), r=nm(r)) %>% filter(q!=r)
pairs <- ava %>% left_join(mash %>% select(q,r,mash_dist=dist), by=c("q","r")) %>%
  mutate(key=ifelse(q<r, paste(q,r), paste(r,q))) %>% group_by(key) %>%
  summarise(a=first(q), b=first(r), ani=mean(ani), mash_dist=mean(mash_dist), .groups="drop")
cat("pairwise ANI: min", round(min(pairs$ani),3), " median", round(median(pairs$ani),3), " max", round(max(pairs$ani),4), "\n")
cat("mash distance: min", signif(min(pairs$mash_dist),3), " median", signif(median(pairs$mash_dist),3), "\n")
# clone criterion: mash distance <= 0.0002 (~99.98% identity) AND ANI >= 99.95
clones <- pairs %>% filter(mash_dist <= 2e-4, ani >= 99.95) %>% arrange(mash_dist)
cat("\ncandidate clone pairs (mash<=2e-4 & ANI>=99.95):", nrow(clones), "\n")
if (nrow(clones)) {
  cl <- clones %>% left_join(ph %>% select(a=accession, sa=strain_name, pa=mic_used, ca=country), by="a") %>%
                   left_join(ph %>% select(b=accession, sb=strain_name, pb=mic_used, cb=country), by="b")
  print(as.data.frame(cl %>% select(a,sa,ca,pa,b,sb,cb,pb,ani,mash_dist)))
}
# transitive clone groups
grp <- setNames(seq_len(nrow(ph)), ph$accession)
for (i in seq_len(nrow(clones))) { g <- min(grp[[clones$a[i]]], grp[[clones$b[i]]]); old <- c(grp[[clones$a[i]]], grp[[clones$b[i]]]); grp[grp %in% old] <- g }
gdf <- data.frame(accession=names(grp), clone_group=as.integer(factor(unlist(grp))), stringsAsFactors=FALSE)
sizes <- gdf %>% count(clone_group, name="group_size")
gdf <- gdf %>% left_join(sizes, by="clone_group")
cat("\nclone groups with >1 member:", sum(sizes$group_size>1), " genomes involved:", sum(sizes$group_size[sizes$group_size>1]), "\n")
# ---- assembly stats
sz <- read.delim(file.path(W,"00_qc/genome_sizes.tsv"), header=FALSE, col.names=c("accession","bases","n_contigs"))
qc <- ph %>% select(accession, strain_name, mic_used, R_call, country, cyp51a_genotype, level, coverage_x, contig_n50) %>%
  left_join(spec, by="accession") %>% left_join(sz, by="accession") %>% left_join(gdf, by="accession") %>%
  mutate(size_Mb=round(bases/1e6,2), species_ok = best_ref=="fumigatus_Af293" & ani_best>=98)
write.csv(qc, file.path(base,"output/table/vrc_phase0_qc.csv"), row.names=FALSE, na="")
cat("\nspecies check: all fumigatus?", all(qc$species_ok), "  ANI-to-Af293 range:", round(min(qc$ani_best),2), "-", round(max(qc$ani_best),2), "\n")
# ---- figures
p1 <- ggplot(pairs, aes(ani)) + geom_histogram(bins=60, fill="grey35") +
  geom_vline(xintercept=99.95, linetype=2, colour="#b2182b") +
  labs(x="Pairwise ANI (%)", y="Genome pairs", title="Pairwise genome identity",
       subtitle="dashed line = clone threshold (99.95%)") + theme_bw(base_size=11)
p2 <- ggplot(qc, aes(size_Mb, contig_n50/1e3, colour=R_call)) + geom_point(size=2.5, alpha=.85) +
  scale_y_log10() + scale_colour_manual(values=c(R="#b2182b", S="#2166ac")) +
  labs(x="Assembly size (Mb)", y="Contig N50 (kb, log)", colour="VRC", title="Assembly quality of the analysis set") + theme_bw(base_size=11)
pdf(file.path(base,"output/figure/vrc_phase0_qc.pdf"), width=11, height=4.5); gridExtra::grid.arrange(p1,p2,ncol=2); dev.off()
cat("\nwrote output/table/vrc_phase0_qc.csv and output/figure/vrc_phase0_qc.pdf\n")
