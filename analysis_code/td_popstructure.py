#!/usr/bin/env python3
"""Detect clonal structure among the 315 A. fumigatus genomes and de-duplicate it.

Constraint estimates assume genomes are independent observations. If a group of isolates is a
recent clonal expansion, genes look conserved in that group simply because there has not been
time to vary - not because selection preserved them. This builds a pairwise distance matrix from
the extracted variants, clusters genomes into lineages, and writes a de-duplicated genome list
(one representative per lineage) for recomputing constraint.
"""
import glob, gzip, os, numpy as np
from collections import defaultdict

TD = "analysis/target_discovery"
files = sorted(glob.glob(os.path.join(TD, "var/*.snp.tsv.gz")))
accs = [os.path.basename(f).replace(".snp.tsv.gz", "") for f in files]
print(f"{len(files)} genomes", flush=True)

# index unique variants seen in >=2 and <=n-2 genomes (informative sites only)
print("indexing variants...", flush=True)
counts = defaultdict(int)
per_genome = []
for f in files:
    s = set()
    with gzip.open(f, "rt") as fh:
        for line in fh:
            c, p, r, a = line.rstrip("\n").split("\t")
            if a in "ACGT" and r in "ACGT":
                k = (c, int(p), a)
                s.add(k); counts[k] += 1
    per_genome.append(s)
n = len(files)
informative = {k for k, v in counts.items() if 2 <= v <= n - 2}
print(f"  {len(counts)} unique alleles, {len(informative)} informative", flush=True)

idx = {k: i for i, k in enumerate(sorted(informative))}
m = len(idx)
M = np.zeros((n, m), dtype=bool)
for gi, s in enumerate(per_genome):
    cols = [idx[k] for k in s if k in idx]
    M[gi, cols] = True
print(f"  matrix {M.shape}", flush=True)

# pairwise Hamming distance (fraction of informative sites that differ)
print("computing pairwise distances...", flush=True)
Mi = M.astype(np.int16)
D = np.zeros((n, n), dtype=np.float32)
for i in range(n):
    diff = np.logical_xor(M[i], M).sum(axis=1)
    D[i] = diff / m
np.fill_diagonal(D, 0.0)

# hierarchical clustering, average linkage
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
Z = linkage(squareform(D, checks=False), method="average")

print("\n=== clonal structure at several thresholds ===")
best = None
for t in [0.02, 0.05, 0.10, 0.15, 0.20, 0.30]:
    cl = fcluster(Z, t=t, criterion="distance")
    sizes = np.bincount(cl)[1:]
    print(f"  cutoff {t:.2f}: {cl.max():3d} lineages, largest {sizes.max():3d}, "
          f"singletons {int((sizes==1).sum()):3d}")
    if t == 0.02:
        best = cl

cl = best
sizes = np.bincount(cl)[1:]
print(f"\nusing cutoff 0.02 -> {cl.max()} lineages")

# representative per lineage: the genome with the most variants (best assembled / least missing)
nvar = np.array([len(s) for s in per_genome])
reps = []
for c in range(1, cl.max() + 1):
    members = np.where(cl == c)[0]
    reps.append(members[np.argmax(nvar[members])])
reps = sorted(reps)

os.makedirs(os.path.join(TD, "popstruct"), exist_ok=True)
with open(os.path.join(TD, "popstruct/lineages.tsv"), "w") as fh:
    fh.write("accession\tlineage\tn_variants\tis_representative\n")
    for i, a in enumerate(accs):
        fh.write(f"{a}\t{cl[i]}\t{nvar[i]}\t{'yes' if i in reps else 'no'}\n")
with open(os.path.join(TD, "popstruct/representatives.txt"), "w") as fh:
    for i in reps:
        fh.write(accs[i] + "\n")
np.save(os.path.join(TD, "popstruct/distance.npy"), D)
print(f"wrote {len(reps)} representatives to popstruct/representatives.txt")

# how clonal is it, really
big = sizes[sizes >= 5]
print(f"\nlineages with >=5 members: {len(big)}, covering {big.sum()} of {n} genomes "
      f"({100*big.sum()/n:.0f}%)")
print(f"largest lineage: {sizes.max()} genomes ({100*sizes.max()/n:.0f}%)")
