#!/usr/bin/env Rscript
# Manhattan/QQ for the genome-wide scans + a consolidated results summary.
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(data.table)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; W <- file.path(base,"analysis/vrc_gwas"); A <- file.path(W,"04_association")
scans <- c(unitigs_bin_dedup46="unitigs, R/S, n=46", unitigs_bin_all55="unitigs, R/S, n=55",
           unitigs_mic_dedup46="unitigs, log2MIC, n=46", snps_bin_dedup46="SNPs, R/S, n=46", snps_mic_dedup46="SNPs, log2MIC, n=46")
summ <- lapply(names(scans), function(s) {
  d <- fread(file.path(A, paste0(s,".tsv")), select=c("variant","af","lrt-pvalue","beta"))
  setnames(d, c("variant","af","p","beta")); d <- d[!is.na(p)]
  np <- length(unique(readLines(file.path(A, paste0(s,"_patterns.txt")))))
  thr <- 0.05/np
  # genomic inflation
  chi <- qchisq(1-d$p, 1); lambda <- median(chi, na.rm=TRUE)/qchisq(0.5,1)
  data.frame(scan=s, label=scans[[s]], tested=nrow(d), patterns=np, threshold=thr,
             n_sig=sum(d$p<thr), n_p1e5=sum(d$p<1e-5), min_p=min(d$p), lambda=round(lambda,3))
}) %>% bind_rows()
write.csv(summ, file.path(base,"output/table/vrc_gwas_scan_summary.csv"), row.names=FALSE)
print(summ, digits=3)
# QQ + Manhattan-style rank plot for the two headline scans
mk <- function(s) {
  d <- fread(file.path(A, paste0(s,".tsv")), select=c("variant","lrt-pvalue")); setnames(d, c("variant","p")); d <- d[!is.na(p)]
  d <- d[order(p)]; n <- nrow(d)
  data.frame(scan=scans[[s]], obs=-log10(d$p[1:min(n,200000)]), exp=-log10(ppoints(n)[1:min(n,200000)]),
             thr=-log10(0.05/length(unique(readLines(file.path(A, paste0(s,"_patterns.txt")))))))
}
qq <- bind_rows(lapply(c("unitigs_bin_dedup46","unitigs_bin_all55","unitigs_mic_dedup46","snps_bin_dedup46"), mk))
p1 <- ggplot(qq, aes(exp, obs)) + geom_abline(slope=1, intercept=0, colour="grey60") +
  geom_point(size=.5, alpha=.5, colour="#2166ac") + geom_hline(aes(yintercept=thr), linetype=2, colour="#b2182b") +
  facet_wrap(~scan, ncol=2) + theme_bw(base_size=11) +
  labs(x=expression(expected~-log[10](p)), y=expression(observed~-log[10](p)),
       title="Voriconazole association scans: no variant reaches genome-wide significance",
       subtitle="red dashed line = Bonferroni threshold on unique presence/absence patterns")
ggsave(file.path(base,"output/figure/vrc_gwas_qq.pdf"), p1, width=9, height=7)
# panel results figure
r <- read.csv(file.path(base,"output/table/vrc_candidate_panel_results.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")) %>%
     filter(test=="substitution", !is.na(p)) %>% mutate(lab=ifelse(p<0.03, variant, NA))
p2 <- ggplot(r, aes(reorder(paste(set,variant), -p), -log10(p), colour=set)) + geom_point(size=1.8, alpha=.85) +
  geom_hline(yintercept=-log10(0.05), linetype=3, colour="grey40") +
  geom_hline(aes(yintercept=-log10(0.05/n_tests)), linetype=2, colour="#b2182b") +
  coord_flip() + facet_wrap(~set, scales="free_y") + theme_bw(base_size=10) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  labs(x="candidate-gene substitutions", y=expression(-log[10](p)), colour=NULL,
       title="Candidate-gene panel: nothing survives multiple-testing correction",
       subtitle="dotted = nominal p 0.05; dashed red = Bonferroni within the test family")
ggsave(file.path(base,"output/figure/vrc_candidate_panel.pdf"), p2, width=10, height=6)
cat("\nwrote vrc_gwas_scan_summary.csv, vrc_gwas_qq.pdf, vrc_candidate_panel.pdf\n")
