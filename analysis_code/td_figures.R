#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table)})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
FIG <- "doc/Aspergillus/target_discovery_2026-08/figures"; dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
d <- fread("doc/Aspergillus/target_discovery_2026-08/tables/constraint_all_genes.csv")
cand <- fread("doc/Aspergillus/target_discovery_2026-08/tables/candidate_targets_ranked.csv")
ACC<-"#33685A"; POS<-"#B2182B"; NEG<-"#2166AC"; GREY<-"#6B7A74"; WARN<-"#8A6410"; RULE<-"#D2DBD6"
op <- function(f,w,h) pdf(file.path(FIG,f), width=w, height=h, useDingbats=FALSE)
ps <- function(...) par(mar=c(4.4,4.6,3.2,1.2), mgp=c(2.6,.7,0), tcl=-.25, cex.axis=.85,
                        cex.lab=.95, las=1, bty="n", ...)
ev <- d[mean_cov_genomes>=250 & tot_seg>=3]

## FIG 1 - the funnel
steps <- c("All protein-coding genes"=9823, "Callable in >=250 genomes"=9445,
           "Evaluable"=9078, "Constrained"=4469, "No human homolog"=2497,
           "No paralog"=1715, "Conserved in >=8 spp"=1240,
           "Not already drugged"=1239, "Length >=100 aa"=1237,
           "Adequately powered"=nrow(cand))
op("fig1_funnel.pdf",7.6,5.0); par(mar=c(4.2,13.5,3.0,3.2), mgp=c(2.5,.6,0), las=1, bty="n", tcl=-.2)
bp <- barplot(rev(steps), horiz=TRUE, col=rev(colorRampPalette(c(ACC,"#9CC4B6"))(length(steps))),
              border=NA, xlim=c(0,10400), cex.names=.78, xlab="genes remaining")
text(rev(steps)+180, bp, rev(steps), adj=0, cex=.75, col=GREY, xpd=NA)
mtext("From 9,823 genes to a ranked shortlist", 3, line=1.5, adj=0, cex=1.0, font=2)
mtext("Each bar is one filter. The largest single cut is human homology - selectivity, not constraint.",
      3, line=.4, adj=0, cex=.72, col=GREY)
dev.off()

## FIG 2 - genome-wide constraint with known targets marked
known <- data.table(
 gene=c("cyp51A","gwt1","DHODH","cyp51B","atrR","hmg1","fks1","hsp90","cnaA","hapE","erg25B","cdr1B"),
 locus=c("AFUA_4G06890","AFUA_1G14870","AFUA_2G11010","AFUA_7G03740","AFUA_2G02690",
         "AFUA_2G03700","AFUA_6G12400","AFUA_5G04170","AFUA_5G09360","AFUA_6G05300",
         "AFUA_4G04820","AFUA_1G14330"))
k <- merge(known, ev, by.x="locus", by.y="locus_tag")
op("fig2_constraint_distribution.pdf",8.0,5.2); ps()
h <- hist(pmin(ev$piNpiS_shrunk,1.5), breaks=60, plot=FALSE)
plot(h, col="#E7ECE8", border="white", main="", xlab=expression(pi[N]/pi[S]~"(shrunk)"),
     ylab="genes", xlim=c(0,1.5))
abline(v=median(ev$piNpiS_shrunk,na.rm=TRUE), col=GREY, lty=2)
text(median(ev$piNpiS_shrunk,na.rm=TRUE), max(h$counts)*.98, " genome median", adj=0, cex=.72, col=GREY)
setorder(k, piNpiS_shrunk)
yy <- seq(max(h$counts)*.86, max(h$counts)*.26, length.out=nrow(k))
for (i in seq_len(nrow(k))) {
  col <- if (k$piNpiS_shrunk[i] < .15) ACC else if (k$piNpiS_shrunk[i] > .5) POS else WARN
  segments(k$piNpiS_shrunk[i], 0, k$piNpiS_shrunk[i], yy[i], col=col, lwd=1.1, lty=3)
  points(k$piNpiS_shrunk[i], yy[i], pch=19, col=col, cex=.75)
  text(k$piNpiS_shrunk[i]+.022, yy[i], k$gene[i], adj=0, cex=.72, col=col, font=2)
}
mtext("Constraint across 9,823 genes, and where the drug targets sit", 3, line=1.5, adj=0, cex=1.0, font=2)
mtext("Green = constrained (resistance-resistant).  Red = permissive - and these are exactly the targets resistance beat.",
      3, line=.4, adj=0, cex=.72, col=GREY)
dev.off()

## FIG 3 - the selectivity paradox
op("fig3_selectivity_paradox.pdf",7.4,5.0); ps()
ev2 <- ev[!is.na(piNpiS_shrunk)][order(piNpiS_shrunk)]
nb <- 20; ev2[, bin := cut(seq_len(.N), nb, labels=FALSE)]
agg <- ev2[, .(frac_human=mean(has_human_homolog), mid=median(piNpiS_shrunk)), by=bin]
plot(agg$bin, 100*agg$frac_human, type="b", pch=19, col=POS, lwd=2, cex=.9,
     xlab="constraint rank bin  (1 = most constrained 5% of genes)",
     ylab="% of genes with a human homolog", ylim=c(20,70), xaxt="n")
axis(1, at=c(1,5,10,15,20), labels=c("1\nmost\nconstrained","5","10","15","20\nleast\nconstrained"),
     padj=.5, cex.axis=.75)
abline(h=100*mean(d$has_human_homolog), col=GREY, lty=2)
text(19.6, 100*mean(d$has_human_homolog)+1.6, "genome baseline 33%", adj=1, cex=.7, col=GREY)
mtext("The most constrained genes are the ones you share with the patient", 3, line=1.5, adj=0, cex=1.0, font=2)
mtext("60% of the 200 most constrained genes have a human homolog, against a 33% genome baseline.",
      3, line=.4, adj=0, cex=.72, col=GREY)
dev.off()

## FIG 4 - top candidates
top <- head(cand, 25)[order(-piNpiS_shrunk)]
op("fig4_top_candidates.pdf",8.4,6.2); par(mar=c(4.4,20,3.2,2.4), mgp=c(2.5,.6,0), las=1, bty="n", tcl=-.2)
lab <- paste0(top$locus_tag, "  ", substr(top$product,1,34))
bp <- barplot(top$piNpiS_shrunk, horiz=TRUE, names.arg=lab, col=ACC, border=NA,
              xlim=c(0,max(top$piNpiS_shrunk)*1.28), cex.names=.62,
              xlab=expression(pi[N]/pi[S]~"  (lower = more constrained)"))
text(top$piNpiS_shrunk+max(top$piNpiS_shrunk)*.02, bp,
     sprintf("%.3f  (%d obs / %.0f exp)", top$piNpiS_shrunk, top$tot_seg, top$exp_tot),
     adj=0, cex=.58, col=GREY, xpd=NA)
mtext("Top 25 candidates after every filter", 3, line=1.5, adj=0, cex=1.0, font=2)
mtext("Constrained · no human homolog · no paralog · pan-Aspergillus · novel",
      3, line=.4, adj=0, cex=.72, col=GREY)
dev.off()

cat("figures written\n"); print(list.files(FIG))

## FIG 5 - constraint vs experimental essentiality (van Rhijn 2024 kinome)
k <- fread("doc/Aspergillus/target_discovery_2026-08/tables/kinome_vanrhijn2024.csv")
k[, locus_tag := toupper(gsub("\\s","",locus_tag))]
mk <- merge(k, d, by="locus_tag")[mean_cov_genomes>=250 & tot_seg>=3]
mk[, essential := homokaryotic_null_generated=="No"]
op("fig5_essentiality_validation.pdf",6.8,5.0); ps()
grp <- list(mk[essential==TRUE]$piNpiS_shrunk, mk[essential==FALSE]$piNpiS_shrunk)
cols2 <- c(ACC, GREY)
plot(NA, xlim=c(.5,2.5), ylim=c(0, quantile(unlist(grp),.97)), xaxt="n",
     xlab="", ylab=expression(pi[N]/pi[S]~"(shrunk)"))
for (i in 1:2) {
  v <- grp[[i]]; set.seed(i)
  points(jitter(rep(i,length(v)),1.5), v, pch=16, cex=.6, col=adjustcolor(cols2[i],.42))
  q <- quantile(v, c(.25,.5,.75))
  rect(i-.24,q[1],i+.24,q[3], border=cols2[i], lwd=1.8)
  segments(i-.24,q[2],i+.24,q[2], col=cols2[i], lwd=3)
  text(i, quantile(unlist(grp),.97)*.98, sprintf("n=%d", length(v)), cex=.78, col=cols2[i])
}
axis(1, at=1:2, labels=c("essential\n(no null obtainable)","non-essential\n(null generated)"),
     padj=.5, cex.axis=.82)
w <- wilcox.test(piNpiS_shrunk ~ essential, data=mk)
r <- rank(-mk$piNpiS_shrunk); n1 <- sum(mk$essential); n0 <- sum(!mk$essential)
auc <- (sum(r[mk$essential]) - n1*(n1+1)/2)/(n1*n0)
mtext("Constraint does predict essentiality - weakly", 3, line=1.6, adj=0, cex=1.0, font=2)
mtext(sprintf("A. fumigatus kinome, van Rhijn 2024.  medians %.3f vs %.3f,  Wilcoxon p = %.4f,  AUC = %.3f",
      median(grp[[1]]), median(grp[[2]]), w$p.value, auc), 3, line=.45, adj=0, cex=.7, col=GREY)
dev.off()
cat("fig5 written\n")
