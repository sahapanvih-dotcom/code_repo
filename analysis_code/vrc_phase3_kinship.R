#!/usr/bin/env Rscript
# Kinship matrix for the LMM: shared root-to-MRCA branch length from the core-SNP ML tree
# (ape::vcv on a phylogeny gives exactly this, i.e. pyseer's phylogeny_distance --lmm matrix).
suppressPackageStartupMessages({library(ape); library(dplyr); library(ggplot2)})
base <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"; W <- file.path(base,"analysis/vrc_gwas")
tr <- read.tree(file.path(W,"03_structure/core_snps.treefile"))
tr <- midpoint <- root(tr, outgroup=tr$tip.label[which.max(diag(vcv(tr)))], resolve.root=TRUE)
K <- vcv(tr); K <- K[order(rownames(K)), order(colnames(K))]
write.table(cbind(data.frame(V1=rownames(K)), as.data.frame(K)), file.path(W,"03_structure/phylogeny_K.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
cat("K:", nrow(K), "x", ncol(K), "\n")
# patristic distance + clustering for lineage assignment
D <- cophenetic(tr); hc <- hclust(as.dist(D), method="average")
for (k in 2:6) assign(paste0("cl",k), cutree(hc, k))
cl <- data.frame(accession=names(cutree(hc,4)), clade_k2=cutree(hc,2), clade_k3=cutree(hc,3), clade_k4=cutree(hc,4))
ph <- read.csv(file.path(base,"output/table/voriconazole_phenotype_fumigatus.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
cy <- read.csv(file.path(base,"output/table/vrc_cyp51A_genotypes.csv"), stringsAsFactors=FALSE, na.strings=c("","NA")) %>%
      select(accession, promoter_genotype, aa_substitutions, clone_group, group_size, known_mechanism)
m <- ph %>% select(accession, strain_name, mic_used, R_call, country, mic_source) %>% left_join(cl, by="accession") %>% left_join(cy, by="accession")
write.csv(m, file.path(base,"output/table/vrc_population_structure.csv"), row.names=FALSE, na="")
cat("\nclade (k=2) vs voriconazole:\n"); print(table(m$clade_k2, m$R_call))
cat("clade (k=4) vs voriconazole:\n"); print(table(m$clade_k4, m$R_call))
cat("clade (k=4) vs country:\n"); print(table(m$clade_k4, m$country))
cat("\nchi-square, clade_k4 vs R/S: p =", signif(suppressWarnings(chisq.test(table(m$clade_k4, m$R_call))$p.value),3), "\n")
# tree figure
pdf(file.path(base,"output/figure/vrc_phylogeny.pdf"), width=9, height=11)
tipcol <- ifelse(m$R_call[match(tr$tip.label, m$accession)]=="R", "#b2182b", "#2166ac")
lab <- paste0(m$strain_name[match(tr$tip.label, m$accession)], "  [", m$mic_used[match(tr$tip.label, m$accession)], "]",
              ifelse(m$promoter_genotype[match(tr$tip.label, m$accession)]!="none", paste0(" ", m$promoter_genotype[match(tr$tip.label, m$accession)]), ""))
tr2 <- tr; tr2$tip.label <- lab
plot(tr2, tip.color=tipcol, cex=.62, no.margin=FALSE, main="A. fumigatus core-SNP phylogeny (n=55)")
legend("bottomleft", legend=c("voriconazole R","voriconazole S"), text.col=c("#b2182b","#2166ac"), bty="n", cex=.8)
mtext("tip label: strain [MIC mg/L] cyp51A promoter genotype", side=1, line=2, cex=.75)
dev.off()
cat("\nwrote phylogeny_K.tsv, output/table/vrc_population_structure.csv, output/figure/vrc_phylogeny.pdf\n")
