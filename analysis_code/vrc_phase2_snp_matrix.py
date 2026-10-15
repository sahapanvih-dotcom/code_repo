#!/usr/bin/env python3
"""Build a biallelic SNP genotype matrix across the 55 A. fumigatus genomes.
Non-R step: interval/genotype assembly over ~10^5 sites x 55 genomes is impractical in base R here.
Genotype coding: 1 = alternate allele, 0 = reference allele (site inside an aligned block),
                 NA = site not covered by any alignment block for that genome.
Outputs a pyseer-compatible presence/absence Rtab plus a core-SNP FASTA alignment for the phylogeny."""
import os, glob, gzip, bisect, sys
from collections import defaultdict

W = "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta/analysis/vrc_gwas"
OUT = f"{W}/02_features"
accs = sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{W}/genomes/*.fna"))

# ---------- coverage intervals ----------
cov = {}
for a in accs:
    d = defaultdict(list)
    for line in open(f"{OUT}/paf/{a}.bed"):
        c, s, e = line.split("\t"); d[c].append((int(s), int(e)))
    merged = {}
    for c, iv in d.items():
        iv.sort(); out = []
        for s, e in iv:
            if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
            else: out.append([s, e])
        merged[c] = (tuple(x[0] for x in out), tuple(x[1] for x in out))
    cov[a] = merged
def covered(a, c, pos0):
    m = cov[a].get(c)
    if not m: return False
    starts, ends = m
    i = bisect.bisect_right(starts, pos0) - 1
    return i >= 0 and pos0 < ends[i]

# ---------- variants ----------
alt = defaultdict(dict)          # (chrom,pos) -> {acc: altbase}
refbase = {}
for a in accs:
    for line in open(f"{OUT}/snps/{a}.vcf"):
        if line.startswith("#"): continue
        f = line.split("\t")
        r, v = f[3], f[4]
        if len(r) != 1 or len(v) != 1 or v == ".": continue      # SNVs only
        key = (f[0], int(f[1]) - 1)
        alt[key][a] = v; refbase[key] = r
print(f"candidate SNV sites: {len(alt)}", flush=True)

# ---------- filter + genotype ----------
MAX_MISS = 0.05
MIN_MAC  = 3
rows, aln = [], {a: [] for a in accs}
kept = 0
for key in sorted(alt):
    c, p = key
    bases = set(alt[key].values())
    if len(bases) != 1: continue                                  # keep strictly biallelic
    b = bases.pop()
    g, miss = [], 0
    for a in accs:
        if a in alt[key]: g.append(1)
        elif covered(a, c, p): g.append(0)
        else: g.append(None); miss += 1
    if miss / len(accs) > MAX_MISS: continue
    mac = min(sum(1 for x in g if x == 1), sum(1 for x in g if x == 0))
    if mac < MIN_MAC: continue
    rows.append((c, p + 1, refbase[key], b, g)); kept += 1
    for a, x in zip(accs, g):
        aln[a].append(b if x == 1 else refbase[key] if x == 0 else "N")
print(f"kept sites (biallelic, <=5% missing, MAC>={MIN_MAC}): {kept}", flush=True)

with gzip.open(f"{OUT}/snp_matrix.Rtab.gz", "wt") as fh:
    fh.write("Variant\t" + "\t".join(accs) + "\n")
    for c, p, r, b, g in rows:
        fh.write(f"{c}_{p}_{r}_{b}\t" + "\t".join("NA" if x is None else str(x) for x in g) + "\n")
with open(f"{OUT}/core_snps.fasta", "w") as fh:
    for a in accs: fh.write(f">{a}\n{''.join(aln[a])}\n")
with open(f"{OUT}/snp_sites.tsv", "w") as fh:
    fh.write("chrom\tpos\tref\talt\tn_alt\n")
    for c, p, r, b, g in rows: fh.write(f"{c}\t{p}\t{r}\t{b}\t{sum(1 for x in g if x==1)}\n")
print("wrote snp_matrix.Rtab.gz, core_snps.fasta, snp_sites.tsv")
