#!/usr/bin/env Rscript
# QC + gene copy-number calling for the C. tropicalis set.
# Copy number = 2 * (gene mean depth / genome-wide median gene depth), since C. tropicalis is diploid.
# Normalising to the median gene rather than the mean makes the estimate robust to a few amplified genes.
# The published ERG11 copy numbers in Fan 2023 are the validation target.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(data.table)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; C <- file.path(base,"analysis/candida")
out <- file.path(base,"doc/Candida")
PANEL <- c(ERG11="CTRG_05283", TAC1="CTRG_05307", UPC2="CTRG_03159", MRR1="CTRG_02269",
           MDR1="CTRG_02412", CDR2="CTRG_02337", ERG3="CTRG_04480", ERG5="CTRG_04983",
           ERG6="CTRG_02396", ERG2="CTRG_04378", ERG24="CTRG_01869", ERG25="CTRG_00908")
fs <- list.files(file.path(C,"depth"), pattern="\\.regions\\.bed\\.gz$", full.names=TRUE)
cat("depth files:", length(fs), "\n")
res <- rbindlist(lapply(fs, function(f) {
  acc <- sub("\\.regions\\.bed\\.gz$","",basename(f))
  d <- fread(f, col.names=c("chrom","start","end","gene","depth"))
  med <- median(d$depth[d$depth > 0])
  g <- d[gene %in% PANEL]
  data.table(accession=acc, median_gene_depth=med,
             gene=names(PANEL)[match(g$gene, PANEL)], locus=g$gene,
             gene_depth=g$depth, copies=round(2*g$depth/med, 2))
}), fill=TRUE)
w <- dcast(res, accession + median_gene_depth ~ gene, value.var="copies")
# QC metrics
qcf <- list.files(file.path(C,"qc"), pattern="\\.fastp\\.json$", full.names=TRUE)
qc <- rbindlist(lapply(qcf, function(f){
  j <- tryCatch(jsonlite::fromJSON(f), error=function(e) NULL); if (is.null(j)) return(NULL)
  data.table(accession=sub("\\.fastp\\.json$","",basename(f)),
             reads_after=j$summary$after_filtering$total_reads,
             q30_rate=round(j$summary$after_filtering$q30_rate,4),
             gc=round(j$summary$after_filtering$gc_content,4),
             dup_rate=round(if(!is.null(j$duplication$rate)) j$duplication$rate else NA,4))
}), fill=TRUE)
m <- merge(w, qc, by="accession", all.x=TRUE)
fwrite(m, file.path(out,"tables/candida_qc_and_copynumber.csv"))
cat("\n=== QC summary (", nrow(m), "isolates so far )\n")
cat("median gene depth: median", round(median(m$median_gene_depth),1),
    " range", round(min(m$median_gene_depth),1), "-", round(max(m$median_gene_depth),1), "\n")
cat("isolates below 10x median depth:", sum(m$median_gene_depth < 10), "\n")
cat("\n=== ERG11 copy number distribution (2 = normal diploid):\n")
print(table(cut(m$ERG11, breaks=c(0,1.5,2.5,3.5,4.5,6.5,100),
                labels=c("<1.5","~2","~3","~4","5-6",">6"))))
cat("\nisolates with ERG11 amplification (>2.5 copies):", sum(m$ERG11 > 2.5, na.rm=TRUE), "\n")
cat("\n=== panel gene copy number, median across isolates:\n")
print(round(sapply(names(PANEL)[names(PANEL) %in% names(m)], function(g) median(m[[g]], na.rm=TRUE)),2))
