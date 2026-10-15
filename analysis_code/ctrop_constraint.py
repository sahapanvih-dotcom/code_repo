#!/usr/bin/env python3
"""Per-gene population-genetic constraint (piN/piS) for *Candida tropicalis*.

Cross-genus replication of the A. fumigatus constraint analysis
(analysis_code/td_constraint.py), computed independently in a second genus from
a completely different data type: short-read VCFs against GCF_000006335.3
(MYA-3404) rather than whole-genome assemblies against Af293.

Differences from the Aspergillus version, all forced by the data:
  - C. tropicalis is DIPLOID, so allele counts are per chromosome: a 0/1 call
    contributes 1 alt allele of 2, a 1/1 call contributes 2 of 2. n = 2 x isolates.
  - Variants come from bcftools VCFs (variant-only), so a position absent from a
    VCF is treated as reference. Callability is therefore taken from the mosdepth
    per-gene depth files rather than from alignment coverage.
  - Clone de-duplication is run here, because this collection is known to contain
    a large clonal Y132F expansion (Fan et al. 2023).

Outputs: analysis/candida/constraint/ctrop_constraint_per_gene.tsv
         analysis/candida/constraint/ctrop_lineage_representatives.txt
"""
import os, gzip, glob, sys, csv
from collections import defaultdict
import numpy as np

BASE = "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
CAN  = os.path.join(BASE, "analysis/candida")
REF  = os.path.join(CAN, "ref/GCF_000006335.3_ASM633v3_genomic.fna")
GFF  = os.path.join(CAN, "ref/genomic.gff")
FAA  = os.path.join(CAN, "ref/protein.faa")
TBL  = os.path.join(BASE, "doc/Candida/tables/ctropicalis_analysis_table.csv")
OUTD = os.path.join(CAN, "constraint")
os.makedirs(OUTD, exist_ok=True)

MIN_QUAL = 30.0
MIN_DP   = 10
DEDUP_SITES = 30000     # sites sampled for the clone-distance matrix
CLONE_CUT   = 0.02      # fraction of sampled sites differing

CODON = {}
# NCBI translation table 12 (alternative yeast nuclear): CTG = Ser, not Leu.
# C. tropicalis is a CUG-clade yeast; with the standard genetic code the
# deposited-protein check below fails for ~63% of genes.
_b, _aa, i = "TCAG", "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", 0
for a in _b:
    for b in _b:
        for c in _b:
            CODON[a + b + c] = _aa[i]; i += 1
COMP = str.maketrans("ACGTN", "TGCAN")
assert CODON["CTG"] == "S" and CODON["CTA"] == "L" and CODON["ATG"] == "M"


def load_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs


# ---------------------------------------------------------------- isolate set
iso = []
with open(TBL) as fh:
    for r in csv.DictReader(fh):
        if r.get("qc_pass", "").upper() == "TRUE":
            iso.append(r["accession"])
iso = [a for a in iso if os.path.exists(os.path.join(CAN, "vcf", a + ".vcf.gz"))]
print(f"{len(iso)} QC-pass isolates with a VCF", flush=True)

# ---------------------------------------------------------------- reference
print("loading reference + GFF...", flush=True)
ref = {k: v.upper() for k, v in load_fasta(REF).items()}
ctg_len = {k: len(v) for k, v in ref.items()}

genes = {}
for line in open(GFF):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "CDS": continue
    at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    pid = at.get("protein_id")
    if not pid: continue
    g = genes.setdefault(pid, {"locus": at.get("locus_tag", ""), "ctg": f[0],
                               "strand": f[6], "exons": [], "product": at.get("product", "")})
    g["exons"].append((int(f[3]) - 1, int(f[4])))
print(f"  {len(genes)} coding transcripts", flush=True)

gene_ids = sorted(genes)
gidx = {p: i for i, p in enumerate(gene_ids)}
gene_arr = {c: np.full(l, -1, dtype=np.int32) for c, l in ctg_len.items()}
cpos_arr = {c: np.full(l, -1, dtype=np.int32) for c, l in ctg_len.items()}
cds_seq = {}
for pid in gene_ids:
    g = genes[pid]
    if g["ctg"] not in ref: continue
    ex = sorted(g["exons"], key=lambda t: t[0], reverse=(g["strand"] == "-"))
    seq, k, c = [], 0, g["ctg"]
    for s, e in ex:
        block = ref[c][s:e]
        idxs = np.arange(e - 1, s - 1, -1) if g["strand"] == "-" else np.arange(s, e)
        if g["strand"] == "-": block = block.translate(COMP)[::-1]
        seq.append(block)
        gene_arr[c][idxs] = gidx[pid]
        cpos_arr[c][idxs] = np.arange(k, k + len(idxs))
        k += len(idxs)
    cds_seq[pid] = "".join(seq)

prot = load_fasta(FAA)
ok = sum(1 for pid, s in cds_seq.items() if pid in prot and
         "".join(CODON.get(s[i:i+3], "X") for i in range(0, len(s) - 2, 3)).rstrip("*") == prot[pid])
print(f"  translation check: {ok}/{len(prot)} transcripts match the deposited protein", flush=True)

# ---------------------------------------------------------------- site counts
site_cache = {}
def syn_sites(codon):
    if codon in site_cache: return site_cache[codon]
    aa = CODON.get(codon)
    if aa is None or aa == "*":
        site_cache[codon] = 0.0; return 0.0
    s = 0.0
    for p in range(3):
        for b in "ACGT":
            if b == codon[p]: continue
            if CODON.get(codon[:p] + b + codon[p+1:]) == aa: s += 1/3
    site_cache[codon] = s
    return s

Nsit, Ssit, ncod = {}, {}, {}
for pid, s in cds_seq.items():
    S, n = 0.0, 0
    for i in range(0, len(s) - 2, 3):
        cod = s[i:i+3]
        if len(cod) < 3 or "N" in cod: continue
        S += syn_sites(cod); n += 1
    Ssit[pid], Nsit[pid], ncod[pid] = S, 3 * n - S, n

# ---------------------------------------------------------------- read VCFs
print("reading VCFs...", flush=True)
per_iso = {}                     # accession -> {(ctg,pos,alt): dosage(1|2)}
for j, acc in enumerate(iso):
    d = {}
    with gzip.open(os.path.join(CAN, "vcf", acc + ".vcf.gz"), "rt") as fh:
        for line in fh:
            if line[0] == "#": continue
            f = line.rstrip("\n").split("\t")
            r, a = f[3], f[4]
            if len(r) != 1 or "," in a or len(a) != 1: continue
            if r not in "ACGT" or a not in "ACGT": continue
            try:
                if float(f[5]) < MIN_QUAL: continue
            except ValueError:
                continue
            fmt = f[8].split(":"); smp = f[9].split(":")
            gt = smp[0].replace("|", "/")
            if gt in ("0/0", "./."): continue
            if "DP" in fmt:
                try:
                    if int(smp[fmt.index("DP")]) < MIN_DP: continue
                except (ValueError, IndexError):
                    pass
            d[(f[0], int(f[1]), a)] = 2 if gt == "1/1" else 1
    per_iso[acc] = d
    if (j + 1) % 50 == 0: print(f"  {j+1}/{len(iso)}", flush=True)

allsites = defaultdict(int)
for d in per_iso.values():
    for key in d: allsites[key] += 1
print(f"  {len(allsites)} unique SNV alleles", flush=True)

# ---------------------------------------------------------------- de-duplicate
print("clone de-duplication...", flush=True)
common = [k for k, c in allsites.items() if 5 <= c <= len(iso) - 5]
rng = np.random.default_rng(0)
if len(common) > DEDUP_SITES:
    common = [common[i] for i in rng.choice(len(common), DEDUP_SITES, replace=False)]
M = np.zeros((len(iso), len(common)), dtype=np.int8)
for i, acc in enumerate(iso):
    d = per_iso[acc]
    for j, k in enumerate(common):
        M[i, j] = d.get(k, 0)
D = np.zeros((len(iso), len(iso)))
for i in range(len(iso)):
    D[i] = (M != M[i]).mean(axis=1)
# average-linkage clustering at CLONE_CUT
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
Z = linkage(squareform((D + D.T) / 2, checks=False), method="average")
cl = fcluster(Z, CLONE_CUT, criterion="distance")
reps, seen = [], set()
for acc, c in zip(iso, cl):
    if c not in seen:
        seen.add(c); reps.append(acc)
print(f"  {len(iso)} isolates -> {len(reps)} lineages at distance {CLONE_CUT}", flush=True)
with open(os.path.join(OUTD, "ctrop_lineage_representatives.txt"), "w") as fh:
    fh.write("\n".join(reps) + "\n")

# ---------------------------------------------------------------- callability
print("per-gene callability from mosdepth...", flush=True)
depth_ok = defaultdict(int)      # locus_tag -> isolates with mean depth >= 10
for acc in iso:
    p = os.path.join(CAN, "depth", acc + ".regions.bed.gz")
    if not os.path.exists(p): continue
    with gzip.open(p, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5: continue
            try:
                if float(f[4]) >= 10: depth_ok[f[3]] += 1
            except ValueError:
                continue

# ---------------------------------------------------------------- classify
def constraint(sample):
    n = 2 * len(sample)
    counts = defaultdict(int)
    for acc in sample:
        for key, dose in per_iso[acc].items():
            counts[key] += dose
    syn_pi, non_pi = defaultdict(float), defaultdict(float)
    syn_seg, non_seg, stop_seg = defaultdict(int), defaultdict(int), defaultdict(int)
    for (c, p, alt), k in counts.items():
        if c not in gene_arr or k == 0: continue
        p0 = p - 1
        if p0 >= len(gene_arr[c]): continue
        gi = gene_arr[c][p0]
        if gi < 0: continue
        pid = gene_ids[gi]
        s = cds_seq.get(pid)
        if s is None: continue
        ci = int(cpos_arr[c][p0]); codon_i, cp = ci // 3, ci % 3
        st = codon_i * 3
        if st + 3 > len(s): continue
        cod = s[st:st+3]
        if "N" in cod: continue
        a = alt if genes[pid]["strand"] == "+" else alt.translate(COMP)
        new = cod[:cp] + a + cod[cp+1:]
        aa0, aa1 = CODON.get(cod), CODON.get(new)
        if aa0 is None or aa1 is None or aa0 == "*": continue
        pi = 2.0 * k * (n - k) / (n * (n - 1)) if n > 1 else 0.0
        if aa0 == aa1:
            syn_seg[pid] += 1; syn_pi[pid] += pi
        else:
            non_seg[pid] += 1; non_pi[pid] += pi
            if aa1 == "*": stop_seg[pid] += 1
    return syn_pi, non_pi, syn_seg, non_seg, stop_seg

print("computing constraint (all isolates)...", flush=True)
A = constraint(iso)
print("computing constraint (lineage representatives)...", flush=True)
B = constraint(reps)

# shrinkage towards the genome-wide mean, as in the Aspergillus analysis
def ratios(res):
    """Return per-gene (piN, piS, counts), the RAW ratio and a shrunk ratio.

    The raw ratio is the primary estimator here. Unlike the A. fumigatus data
    (median 14 segregating sites per gene) this collection carries a median of
    ~108 sites per gene, so the raw ratio is already stable; the genome-mean
    shrinkage that stabilises the Aspergillus estimates instead compresses these
    ones towards the genome median and erases real signal (with shrinkage,
    ribosomal proteins become indistinguishable from the genome median; with the
    raw ratio they are 2x more constrained, as in A. fumigatus).
    """
    syn_pi, non_pi, syn_seg, non_seg, stop_seg = res
    out = {}
    for pid in cds_seq:
        ps = syn_pi[pid] / Ssit[pid] if Ssit[pid] > 0 else 0.0
        pn = non_pi[pid] / Nsit[pid] if Nsit[pid] > 0 else 0.0
        out[pid] = (pn, ps, syn_seg[pid], non_seg[pid], stop_seg[pid])
    gpn = np.mean([v[0] for v in out.values()])
    gps = np.mean([v[1] for v in out.values()])
    raw, sh = {}, {}
    for pid, (pn, ps, ss, ns, st) in out.items():
        raw[pid] = (pn / ps) if ps > 0 else np.nan
        sh[pid] = ((pn + gpn) / (ps + gps)) if (ps + gps) > 0 else np.nan
    return out, sh, raw

rawA, shA, ratA = ratios(A)
rawB, shB, ratB = ratios(B)

out = os.path.join(OUTD, "ctrop_constraint_per_gene.tsv")
with open(out, "w") as fh:
    fh.write("protein_id\tlocus_tag\tproduct\tprotein_len\tsyn_sites\tnonsyn_sites\t"
             "syn_seg\tnonsyn_seg\tstop_seg\tpi_syn\tpi_nonsyn\tpiNpiS\tpiNpiS_dedup\t"
             "piNpiS_shrunk\tpiNpiS_shrunk_dedup\tn_iso_callable\n")
    for pid in gene_ids:
        if pid not in cds_seq: continue
        g = genes[pid]
        pn, ps, ss, ns, st = rawA[pid]
        fh.write(f"{pid}\t{g['locus']}\t{g['product']}\t{len(cds_seq[pid])//3}\t"
                 f"{Ssit[pid]:.1f}\t{Nsit[pid]:.1f}\t{ss}\t{ns}\t{st}\t"
                 f"{ps:.6g}\t{pn:.6g}\t{ratA[pid]:.6g}\t{ratB[pid]:.6g}\t"
                 f"{shA[pid]:.6g}\t{shB[pid]:.6g}\t"
                 f"{depth_ok.get(g['locus'], 0)}\n")
print("wrote", out, flush=True)
print(f"isolates {len(iso)}; lineages {len(reps)}; "
      f"median raw piN/piS {np.nanmedian(list(ratA.values())):.3f}", flush=True)
