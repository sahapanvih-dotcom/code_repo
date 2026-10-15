#!/usr/bin/env Rscript
# Classify NCBI Aspergillus assemblies as suitable (PASS) / unsuitable (FAIL) for a
# MIC-association (GWAS-style variant) analysis, then relocate accession folders.
# PASS requires ALL of:
#   depth >= 30x  (or depth unreported AND contig N50 >= 100 kb: classic Sanger/long-read references)
#   total length within +/-15 % of species expected size
#   contig N50 >= 50 kb;  <= 2000 contigs;  not flagged atypical by NCBI
#   not a GenBank (GCA) duplicate of a RefSeq (GCF) copy already in the set
suppressPackageStartupMessages({library(ggplot2)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
m <- read.csv(file.path(base,"output/table/ncbi_assembly_metrics.csv"), stringsAsFactors=FALSE)
m$coverage_x <- suppressWarnings(as.numeric(m$coverage_x))
m$size_ratio <- m$total_len_Mb / m$expected_Mb
gcf_ids <- sub("^GCF_","",m$accession[grepl("^GCF_",m$accession)])
m$dup_of_refseq <- grepl("^GCA_",m$accession) & sub("^GCA_","",m$accession) %in% gcf_ids

reasons <- with(m, {
  r <- rep("", nrow(m))
  add <- function(r, cond, txt) { cond[is.na(cond)] <- FALSE; ifelse(cond, ifelse(r=="", txt, paste(r, txt, sep="; ")), r) }
  r <- add(r, !is.na(coverage_x) & coverage_x < 30, "depth<30x")
  r <- add(r, is.na(coverage_x) & contig_n50 < 1e5, "depth not reported")
  r <- add(r, abs(size_ratio-1) > 0.15, "size outside +/-15%")
  r <- add(r, contig_n50 < 5e4, "contigN50<50kb")
  r <- add(r, n_contigs > 2000, ">2000 contigs")
  r <- add(r, atypical=="True", "NCBI atypical")
  r <- add(r, dup_of_refseq, "GCA duplicate of RefSeq copy")
  r
})
m$status <- ifelse(reasons=="", "PASS", "FAIL"); m$fail_reason <- reasons
m$dest <- ifelse(m$status=="PASS","coverage_pass","coverage_fail")
out <- m[order(m$species, m$status, -ifelse(is.na(m$coverage_x),0,m$coverage_x)), ]
write.csv(out, file.path(base,"output/table/assembly_coverage_filter.csv"), row.names=FALSE)
summ <- as.data.frame.matrix(table(m$species, m$status)); summ$species <- rownames(summ); summ <- summ[,c("species","PASS","FAIL")]
write.csv(summ, file.path(base,"output/table/assembly_coverage_filter_summary.csv"), row.names=FALSE)
print(summ); cat("TOTAL PASS", sum(m$status=="PASS"), "FAIL", sum(m$status=="FAIL"), "\n")
print(sort(table(unlist(strsplit(m$fail_reason[m$status=="FAIL"],"; "))), decreasing=TRUE))

# Figure: depth vs contig N50 per species
pd <- m; pd$coverage_plot <- ifelse(is.na(pd$coverage_x), 1, pmax(pd$coverage_x,1)); pd$contig_n50 <- pmax(pd$contig_n50,100)
p <- ggplot(pd, aes(coverage_plot, contig_n50, colour=status)) +
  geom_hline(yintercept=5e4, linetype=2, colour="grey50") + geom_vline(xintercept=30, linetype=2, colour="grey50") +
  geom_point(alpha=.6, size=1.4) + scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=c(PASS="#1b7837", FAIL="#b2182b")) +
  facet_wrap(~species, ncol=4) + theme_bw(base_size=11) +
  labs(x="Reported sequencing depth (x; unreported plotted at 1)", y="Contig N50 (bp)", colour=NULL)
ggsave(file.path(base,"output/figure/assembly_coverage_filter.pdf"), p, width=12, height=8)

# Relocate accession folders
for (i in seq_len(nrow(m))) {
  sp <- paste0("Aspergillus_", m$species[i]); src <- file.path(base,"raw_data",sp,m$accession[i])
  dstdir <- file.path(base,"raw_data",m$dest[i],sp); dir.create(dstdir, recursive=TRUE, showWarnings=FALSE)
  if (dir.exists(src)) stopifnot(file.rename(src, file.path(dstdir, m$accession[i])))
}
# species-level NCBI report files -> raw_data/ncbi_reports/<species>/
for (sp in unique(paste0("Aspergillus_", m$species))) {
  d <- file.path(base,"raw_data",sp); rd <- file.path(base,"raw_data","ncbi_reports",sp); dir.create(rd, recursive=TRUE, showWarnings=FALSE)
  for (f in list.files(d, full.names=TRUE)) file.rename(f, file.path(rd, basename(f)))
  if (length(list.files(d, all.files=TRUE, no..=TRUE))==0) unlink(d, recursive=TRUE)
}
cat("relocation done\n")
