#!/usr/bin/env Rscript
# Figures for the Candida C. tropicalis fluconazole report. Output: PDF (house convention).
suppressPackageStartupMessages({library(data.table)})
B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; setwd(B)
FIG <- file.path(B,"doc/Candida/figures")
a  <- fread("doc/Candida/tables/ctropicalis_analysis_table.csv")[qc_pass==TRUE]
ph <- fread("doc/Candida/tables/ctropicalis_fluconazole_phenotype.csv")
ACC <- "#33685A"; POS <- "#B2182B"; NEG <- "#2166AC"; GREY <- "#6B7A74"; RULE <- "#D2DBD6"
pdf_open <- function(f,w,h) pdf(file.path(FIG,f), width=w, height=h, useDingbats=FALSE)
par_std <- function(...) par(mar=c(4.2,4.4,2.4,1.2), mgp=c(2.5,.7,0), tcl=-.25,
                             cex.axis=.85, cex.lab=.95, las=1, bty="n", ...)

## FIG 1 - validation against Fan 2023 published copy number
m <- merge(a, ph[,.(accession=strain_name, resistance_genotype)], by="accession")
m <- m[grepl("copy number", resistance_genotype)]
m[, pub_cn := as.numeric(sub(".*copy number ([0-9]+).*","\\1",resistance_genotype))]
m <- m[!is.na(pub_cn)]
pdf_open("fig1_cnv_validation.pdf",6.6,5.2); par_std()
set.seed(1)
plot(jitter(m$pub_cn,.6), m$CN_ERG11, pch=21, cex=.85, lwd=.7,
     bg=ifelse(m$R_call=="R", adjustcolor(POS,.55), adjustcolor(NEG,.45)), col="white",
     xlab="Published ERG11 copy number (Fan et al. 2023)",
     ylab="Copy number called from read depth", xlim=c(1,14), ylim=c(1,14))
abline(0,1,col=GREY,lty=2,lwd=1.1); abline(h=2.5,col=ACC,lty=3,lwd=1)
med <- m[, .(y=median(CN_ERG11)), by=pub_cn][order(pub_cn)]
points(med$pub_cn, med$y, pch=18, cex=1.5, col=ACC); lines(med$pub_cn, med$y, col=ACC, lwd=1.6)
legend("topleft", bty="n", cex=.8,
       legend=c("resistant isolate","susceptible isolate","median of our calls","1:1 agreement","our amplification threshold (2.5)"),
       pch=c(21,21,18,NA,NA), lty=c(NA,NA,1,2,3), lwd=c(NA,NA,1.6,1.1,1),
       pt.bg=c(adjustcolor(POS,.55),adjustcolor(NEG,.45),NA,NA,NA), col=c("white","white",ACC,GREY,ACC))
mtext(sprintf("n = %d isolates   Pearson r = %.3f   sensitivity 100%%   specificity 91%%",
      nrow(m), cor(m$CN_ERG11,m$pub_cn)), 3, line=.6, adj=0, cex=.8, col=GREY)
dev.off()

## FIG 2 - MIC by mechanism class
lv <- c("neither","amplification only","Y132F only","Y132F + amplification")
a[, mech := factor(mechanism, levels=lv)]
cols <- c(GREY, NEG, "#8A6410", POS)
pdf_open("fig2_mechanism_classes.pdf",7.4,5.2); par_std(mar=c(6.4,4.4,2.4,1.2))
bx <- boxplot(log2_mic ~ mech, data=a, plot=FALSE)
plot(NA, xlim=c(.5,4.5), ylim=range(a$log2_mic)+c(-.6,.9), xaxt="n", yaxt="n",
     xlab="", ylab="Fluconazole MIC (mg/L)")
abline(h=log2(8), col=POS, lty=3, lwd=1)
for (i in 1:4) {
  d <- a[mech==lv[i]]
  set.seed(i); points(jitter(rep(i,nrow(d)),1.6), d$log2_mic, pch=16, cex=.5,
                      col=adjustcolor(cols[i],.35))
  q <- quantile(d$log2_mic, c(.25,.5,.75))
  rect(i-.26,q[1],i+.26,q[3], border=cols[i], lwd=1.6, col=NA)
  segments(i-.26,q[2],i+.26,q[2], col=cols[i], lwd=2.8)
  text(i, max(a$log2_mic)+.85, sprintf("n=%d\n%d R", nrow(d), sum(d$R_call=="R")), cex=.72, col=cols[i])
}
tk <- c(-3,-1,1,3,5,7,9); axis(2, at=tk, labels=2^tk)
axis(1, at=1:4, labels=FALSE)
text(1:4, par("usr")[3]-.9, srt=25, adj=1, xpd=NA, cex=.82,
     labels=c("neither","amplification\nonly","Y132F\nonly","Y132F +\namplification"))
text(4.45, log2(8), "R breakpoint", col=POS, cex=.68, adj=1, pos=3, offset=.2)
mtext("Amplification alone shifts nothing; on a Y132F background it adds ~7 dilutions",
      3, line=.6, adj=0, cex=.8, col=GREY)
dev.off()

## FIG 3 - the interaction, stratified
pdf_open("fig3_interaction_stratified.pdf",8.2,4.6)
par(mfrow=c(1,2), mar=c(4.3,4.3,3.2,1), mgp=c(2.5,.7,0), tcl=-.25, cex.axis=.82, cex.lab=.9, las=1, bty="n")
for (y in c(0,1)) {
  d <- a[Y132F==y]
  plot(d$CN_ERG11, d$log2_mic, pch=21, cex=.8, lwd=.6, col="white",
       bg=ifelse(d$R_call=="R", adjustcolor(POS,.55), adjustcolor(NEG,.45)),
       xlab="ERG11 copy number", ylab="Fluconazole MIC (mg/L)",
       xlim=c(1,14), ylim=range(a$log2_mic), yaxt="n")
  axis(2, at=tk, labels=2^tk)
  abline(v=2.5, col=ACC, lty=3); abline(h=log2(8), col=POS, lty=3)
  fit <- lm(log2_mic ~ CN_ERG11, data=d); sl_p <- summary(fit)$coefficients[2,4]
  xx <- seq(min(d$CN_ERG11), max(d$CN_ERG11), length=50)
  # dashed grey when the trend is not significant, so the drawn line cannot mislead
  lines(xx, predict(fit, data.frame(CN_ERG11=xx)),
        col=if (sl_p<.05) ACC else GREY, lwd=if (sl_p<.05) 2 else 1.3, lty=if (sl_p<.05) 1 else 2)
  st <- if (y==0) "Y132F absent" else "Y132F present"
  bt <- round(coef(lm(log2_mic ~ ERG11_amplified, data=d))[2],2)
  pv <- format(wilcox.test(log2_mic ~ ERG11_amplified, data=d)$p.value, digits=2)
  rho <- round(cor(d$CN_ERG11, d$log2_mic, method="spearman"),2)
  mtext(sprintf("%s   (n = %d)", st, nrow(d)), 3, line=1.9, adj=0, cex=.92, font=2)
  mtext(sprintf("amplified vs not:  beta = %+.2f log2,  p = %s", bt, pv),
        3, line=.85, adj=0, cex=.72, col=if (y==0) GREY else POS)
  mtext(sprintf("copy number vs MIC:  Spearman rho = %.2f,  slope p = %s", rho, format(sl_p, digits=2)),
        3, line=.1, adj=0, cex=.72, col=GREY)
}
dev.off()

## FIG 4 - cross-genus: the two shared residues
asp <- fread("output/table/vrc_cyp51A_genotypes.csv")[!is.na(mic_used)]
asp[, log2mic := log2(mic_used)]
asp[, Y121F := grepl("Y121F", aa_substitutions)]
asp[, G448S := grepl("G448S", aa_substitutions)]
pdf_open("fig4_shared_residues.pdf",8.2,5.1)
par(mfrow=c(1,2), mar=c(4.6,4.4,3.4,1), oma=c(1.4,0,0,0), mgp=c(2.5,.7,0), tcl=-.25, cex.axis=.82, cex.lab=.9, las=1, bty="n")
panel <- function(vals, groups, labs, ttl, sub, ylab) {
  ymin <- min(vals)-.5; ymax <- max(vals)+1.4
  plot(NA, xlim=c(.5,length(labs)+.5), ylim=c(ymin,ymax), xaxt="n", yaxt="n", xlab="", ylab=ylab)
  tk2 <- seq(floor(ymin), ceiling(ymax), by=2); axis(2, at=tk2, labels=2^tk2)
  for (i in seq_along(labs)) {
    v <- vals[groups==i]
    if (!length(v)) next
    set.seed(i); points(jitter(rep(i,length(v)),1.4), v, pch=16, cex=.55,
                        col=adjustcolor(if (i%%2==0) POS else GREY,.4))
    q <- quantile(v, c(.25,.5,.75))
    rect(i-.24,q[1],i+.24,q[3], border=if (i%%2==0) POS else GREY, lwd=1.5)
    segments(i-.24,q[2],i+.24,q[2], col=if (i%%2==0) POS else GREY, lwd=2.6)
    text(i, ymax-.2, sprintf("n=%d", length(v)), cex=.7, col=GREY)
  }
  axis(1, at=seq_along(labs), labels=labs, cex.axis=.75)
  mtext(ttl, 3, line=1.6, adj=0, cex=.92, font=2)
  mtext(sub, 3, line=.4, adj=0, cex=.75, col=GREY)
}
v1 <- c(a$log2_mic, asp$log2mic)
g1 <- c(ifelse(a$Y132F==1,2,1), ifelse(asp$Y121F,4,3))
panel(v1, g1, c("Ct\nwild-type","Ct\nY132F","Af\nwild-type","Af\nY121F"),
      "Shared residue 1  -  Ct Tyr132 = Af Tyr121",
      "same tyrosine, same F substitution, both genera", "MIC (mg/L)")
v2 <- c(a$log2_mic, asp$log2mic)
g2 <- c(ifelse(a$G464S==1,2,1), ifelse(asp$G448S,4,3))
panel(v2, g2, c("Ct\nwild-type","Ct\nG464S","Af\nwild-type","Af\nG448S"),
      "Shared residue 2  -  Ct Gly464 = Af Gly448",
      "same glycine, same S substitution, both genera", "MIC (mg/L)")
mtext("Ct = C. tropicalis, fluconazole MIC (n=424).  Af = A. fumigatus, voriconazole MIC (n=55).  Different drugs - compare the within-genus contrast, not the absolute values.",
      1, line=2.6, adj=1, outer=TRUE, cex=.62, col=GREY)
dev.off()

## FIG 5 - residue conservation map along the protein
rm_ <- fread("doc/Candida/tables/erg11_cyp51A_residue_map.csv")
pdf_open("fig5_residue_map.pdf",9.0,4.0); par(mar=c(3.8,5.6,3.4,1.4), bty="n")
plot(NA, xlim=c(0,530), ylim=c(0,3), xaxt="n", yaxt="n", xlab="", ylab="")
rect(1,1.95,528,2.40, col="#E7ECE8", border=RULE)
text(-14,2.17,"Ct ERG11", adj=1, xpd=NA, cex=.82, font=2)
text(-14,1.96,"528 aa", adj=1, xpd=NA, cex=.66, col=GREY)
rect(1,0.70,515,1.15, col="#E7ECE8", border=RULE)
text(-14,0.92,"Af cyp51A", adj=1, xpd=NA, cex=.82, font=2)
text(-14,0.71,"515 aa", adj=1, xpd=NA, cex=.66, col=GREY)
ct <- rm_[direction=="Ct->Af"]; af <- rm_[direction=="Af->Ct"]
shared_ct <- c("Y132F","G464S"); shared_af <- c("Y121F","G448S")
for (i in 1:nrow(ct)) {
  p <- as.numeric(gsub("\\D","",ct$candida_res[i]))
  q <- suppressWarnings(as.numeric(gsub("\\D","",ct$asp_res[i])))
  sh <- ct$candida_allele[i] %in% shared_ct
  col <- if (sh) POS else if (ct$same_residue[i]) ACC else GREY
  segments(p,1.95,ifelse(is.na(q),p,q),1.15, col=adjustcolor(col, if (sh) .9 else .3), lwd=if (sh) 2.4 else .9)
  points(p,2.40, pch=25, bg=col, col=col, cex=if (sh) 1.0 else .6, xpd=NA)
}
for (i in 1:nrow(af)) {
  q <- as.numeric(gsub("\\D","",af$asp_res[i]))
  sh <- af$candida_allele[i] %in% shared_af
  points(q,0.70, pch=24, bg=if (sh) POS else NEG, col=if (sh) POS else NEG, cex=if (sh) 1.0 else .6, xpd=NA)
}
# shared positions: one clean label per protein, above/below
for (i in which(ct$candida_allele %in% shared_ct)) {
  p <- as.numeric(gsub("\\D","",ct$candida_res[i]))
  text(p,2.66, ct$candida_allele[i], cex=.78, col=POS, font=2, xpd=NA)
}
for (i in which(af$candida_allele %in% shared_af)) {
  q <- as.numeric(gsub("\\D","",af$asp_res[i]))
  text(q,0.40, af$candida_allele[i], cex=.78, col=POS, font=2, xpd=NA)
}
# Aspergillus-only alleles: staggered two rows so labels never collide
afo <- af[!(candida_allele %in% shared_af)]
afo[, pos := as.numeric(gsub("\\D","",asp_res))]
setorder(afo, pos)
shared_pos <- as.numeric(gsub("\\D","", af[candida_allele %in% shared_af]$asp_res))
for (i in seq_len(nrow(afo))) {
  yy <- if (i %% 2 == 1) 0.40 else 0.16
  # never sit on the row used by a nearby shared-residue label
  if (yy == 0.40 && min(abs(afo$pos[i] - shared_pos)) < 32) yy <- 0.16
  segments(afo$pos[i], 0.62, afo$pos[i], yy+0.07, col=adjustcolor(NEG,.5), lwd=.7)
  text(afo$pos[i], yy, afo$candida_allele[i], cex=.66, col=NEG, xpd=NA)
}
axis(1, at=seq(0,500,100), cex.axis=.75, line=-0.6)
mtext("residue position", 1, line=1.9, cex=.8, col=GREY)
mtext("Where each genus mutates a shared drug target", 3, line=1.7, adj=0, cex=.98, font=2)
mtext("49.4% identity over 508 aligned residues.  Red = mutated in BOTH genera.  Green = residue conserved, only Candida uses it.  Blue = only Aspergillus.",
      3, line=.5, adj=0, cex=.66, col=GREY)
dev.off()

cat("figures written to", FIG, "\n"); print(list.files(FIG))
