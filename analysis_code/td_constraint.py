#!/usr/bin/env python3
"""Per-gene population-genetic constraint for A. fumigatus from 315 clinical/environmental genomes.

For every Af293 protein-coding gene this computes, from assembly-to-reference alignments:
  - piN/piS   nucleotide diversity at nonsynonymous vs synonymous sites (Nei-Gojobori site counting).
              Within-species polymorphism, so piN/piS - NOT dN/dS, which is between-species divergence.
  - pN/pS     the same ratio computed on counts of segregating sites (less robust; reported for comparison)
  - nonsense  stop-gained variants, a direct loss-of-function-intolerance signal
  - callability: how many of the 315 genomes actually aligned across the gene. Absence of variants in
              an unaligned gene is absence of DATA, not evidence of constraint, so genes with poor
              coverage are flagged and excluded downstream.
"""
import os, gzip, glob, sys
from collections import defaultdict
import numpy as np

BASE = "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
TD   = os.path.join(BASE, "analysis/target_discovery")
REF  = os.path.join(TD, "ref/Af293.fna")
GFF  = os.path.join(TD, "ref/Af293.gff")
OUT  = os.path.join(TD, "constraint_per_gene.tsv")
# optional: restrict to a de-duplicated genome subset (population-structure correction)
SUBSET = os.environ.get("TD_SUBSET")           # path to a file of accessions, one per line
if SUBSET:
    OUT = os.environ.get("TD_OUT", OUT)

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in _b:
    for b in _b:
        for c in _b:
            CODON[a+b+c] = _aa[i]; i += 1
COMP = str.maketrans("ACGTN", "TGCAN")

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

print("loading reference...", flush=True)
ref = {k: v.upper() for k, v in load_fasta(REF).items()}
ctg_len = {k: len(v) for k, v in ref.items()}

print("parsing GFF...", flush=True)
genes = {}   # protein_id -> dict
for line in open(GFF):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "CDS": continue
    at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    pid = at.get("protein_id")
    if not pid: continue
    g = genes.setdefault(pid, {"locus": at.get("locus_tag", ""), "ctg": f[0],
                               "strand": f[6], "exons": [],
                               "product": at.get("product", "")})
    g["exons"].append((int(f[3]) - 1, int(f[4])))     # 0-based half-open
print(f"  {len(genes)} protein-coding transcripts", flush=True)

# ---- build coordinate maps and reference CDS -------------------------------
gene_ids = sorted(genes)
gidx = {p: i for i, p in enumerate(gene_ids)}
gene_arr = {c: np.full(l, -1, dtype=np.int32) for c, l in ctg_len.items()}
cpos_arr = {c: np.full(l, -1, dtype=np.int32) for c, l in ctg_len.items()}
cds_seq, bad = {}, 0
for pid in gene_ids:
    g = genes[pid]
    ex = sorted(g["exons"], key=lambda t: t[0], reverse=(g["strand"] == "-"))
    seq, k, c = [], 0, g["ctg"]
    if c not in ref: bad += 1; continue
    for s, e in ex:
        block = ref[c][s:e]
        if g["strand"] == "-":
            block = block.translate(COMP)[::-1]
            idxs = np.arange(e - 1, s - 1, -1)
        else:
            idxs = np.arange(s, e)
        seq.append(block)
        gene_arr[c][idxs] = gidx[pid]
        cpos_arr[c][idxs] = np.arange(k, k + len(idxs))
        k += len(idxs)
    cds_seq[pid] = "".join(seq)
print(f"  built CDS for {len(cds_seq)} transcripts ({bad} skipped)", flush=True)

# sanity: translation should match the deposited proteins
prot = load_fasta(os.path.join(TD, "ref/Af293.faa"))
ok = 0
for pid, s in list(cds_seq.items()):
    if pid not in prot: continue
    tr = "".join(CODON.get(s[i:i+3], "X") for i in range(0, len(s) - 2, 3)).rstrip("*")
    if tr == prot[pid]: ok += 1
print(f"  translation check: {ok}/{len(prot)} transcripts match the deposited protein", flush=True)

# ---- Nei-Gojobori site counting -------------------------------------------
print("counting synonymous / nonsynonymous sites...", flush=True)
site_cache = {}
def syn_sites(codon):
    if codon in site_cache: return site_cache[codon]
    aa = CODON.get(codon)
    if aa is None or aa == "*": site_cache[codon] = 0.0; return 0.0
    s = 0.0
    for p in range(3):
        for b in "ACGT":
            if b == codon[p]: continue
            alt = codon[:p] + b + codon[p+1:]
            if CODON.get(alt) == aa: s += 1/3
    site_cache[codon] = s
    return s

Nsit, Ssit, ncod = {}, {}, {}
for pid, s in cds_seq.items():
    S = 0.0; n = 0
    for i in range(0, len(s) - 2, 3):
        cod = s[i:i+3]
        if len(cod) < 3 or "N" in cod: continue
        S += syn_sites(cod); n += 1
    Ssit[pid] = S; Nsit[pid] = 3*n - S; ncod[pid] = n

# ---- aggregate variants across genomes ------------------------------------
print("aggregating variants across genomes...", flush=True)
snp_files = sorted(glob.glob(os.path.join(TD, "var/*.snp.tsv.gz")))
if SUBSET:
    keep = {l.strip() for l in open(SUBSET) if l.strip()}
    snp_files = [f for f in snp_files
                 if os.path.basename(f).replace(".snp.tsv.gz","") in keep]
    print(f"  SUBSET active: {len(snp_files)} genomes", flush=True)
ngen = len(snp_files)
counts = defaultdict(int)
for j, f in enumerate(snp_files):
    with gzip.open(f, "rt") as fh:
        for line in fh:
            c, p, r, a = line.rstrip("\n").split("\t")
            if a not in "ACGT" or r not in "ACGT": continue
            counts[(c, int(p), a)] += 1
    if (j+1) % 50 == 0: print(f"  {j+1}/{ngen}", flush=True)
print(f"  {len(counts)} unique variant alleles across {ngen} genomes", flush=True)

# ---- coverage per position ------------------------------------------------
print("computing callability...", flush=True)
cov = {c: np.zeros(l, dtype=np.int16) for c, l in ctg_len.items()}
cov_files = sorted(glob.glob(os.path.join(TD, "var/*.cov.bed")))
if SUBSET:
    cov_files = [f for f in cov_files
                 if os.path.basename(f).replace(".cov.bed","") in keep]
for f in cov_files:
    seen = {c: np.zeros(l, dtype=bool) for c, l in ctg_len.items()}
    for line in open(f):
        c, s, e = line.split()
        if c in seen: seen[c][int(s):int(e)] = True
    for c in cov: cov[c] += seen[c]

# ---- classify --------------------------------------------------------------
print("classifying variants...", flush=True)
syn_pi   = defaultdict(float); non_pi  = defaultdict(float)
syn_seg  = defaultdict(int);   non_seg = defaultdict(int); stop_seg = defaultdict(int)
n = ngen
for (c, p, alt), k in counts.items():
    if c not in gene_arr: continue
    gi = gene_arr[c][p]
    if gi < 0: continue
    pid = gene_ids[gi]
    s = cds_seq.get(pid)
    if s is None: continue
    ci = int(cpos_arr[c][p])
    codon_i, cp = ci // 3, ci % 3
    st = codon_i * 3
    if st + 3 > len(s): continue
    cod = s[st:st+3]
    if "N" in cod: continue
    a = alt if genes[pid]["strand"] == "+" else alt.translate(COMP)
    new = cod[:cp] + a + cod[cp+1:]
    aa0, aa1 = CODON.get(cod), CODON.get(new)
    if aa0 is None or aa1 is None or aa0 == "*": continue
    # nucleotide diversity contribution for a haploid sample of n genomes
    pi = 2.0 * k * (n - k) / (n * (n - 1)) if n > 1 else 0.0
    if aa0 == aa1:
        syn_seg[pid] += 1; syn_pi[pid] += pi
    else:
        non_seg[pid] += 1; non_pi[pid] += pi
        if aa1 == "*": stop_seg[pid] += 1

# ---- write -----------------------------------------------------------------
print("writing table...", flush=True)
with open(OUT, "w") as fh:
    fh.write("protein_id\tlocus_tag\tcontig\tstrand\tn_codons\tprotein_len\t"
             "syn_sites\tnonsyn_sites\tsyn_seg\tnonsyn_seg\tstop_seg\t"
             "pi_syn\tpi_nonsyn\tpiN_piS\tpN_pS\tmean_cov_genomes\tmin_cov_genomes\tproduct\n")
    for pid in gene_ids:
        if pid not in cds_seq: continue
        g = genes[pid]; c = g["ctg"]
        idxs = np.where(gene_arr[c] == gidx[pid])[0]
        mc = float(cov[c][idxs].mean()) if len(idxs) else 0.0
        mn = int(cov[c][idxs].min()) if len(idxs) else 0
        S, N = Ssit[pid], Nsit[pid]
        piS = syn_pi[pid] / S if S > 0 else float("nan")
        piN = non_pi[pid] / N if N > 0 else float("nan")
        ratio = piN / piS if piS and piS > 0 else float("nan")
        pS = syn_seg[pid] / S if S > 0 else float("nan")
        pN = non_seg[pid] / N if N > 0 else float("nan")
        pr = pN / pS if pS and pS > 0 else float("nan")
        fh.write(f"{pid}\t{g['locus']}\t{c}\t{g['strand']}\t{ncod[pid]}\t{len(cds_seq[pid])//3}\t"
                 f"{S:.1f}\t{N:.1f}\t{syn_seg[pid]}\t{non_seg[pid]}\t{stop_seg[pid]}\t"
                 f"{piS:.6g}\t{piN:.6g}\t{ratio:.4f}\t{pr:.4f}\t{mc:.1f}\t{mn}\t{g['product']}\n")
print("wrote", OUT, flush=True)
