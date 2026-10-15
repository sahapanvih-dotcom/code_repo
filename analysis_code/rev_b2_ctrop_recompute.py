#!/usr/bin/env python3
"""rev_b2: recompute C. tropicalis per-gene piN/piS under matched settings and
under isolate subsampling, reusing the project's existing implementation
(analysis_code/ctrop_constraint.py) verbatim: same translation table 12, same
GFF/CDS build, same Nei-Gojobori site counting, same diploid dosage handling,
same pi = 2k(n-k)/(n(n-1)) with n = 2 x isolates, same VCF quality gates
(QUAL >= 30, DP >= 10).

Reason for Python: it reuses the project's variant parser, as instructed.
All statistics/plots are in R.

Stages:
  cache : parse all QC-pass VCFs once -> per-isolate sparse dosage index (.npz)
  agg   : aggregate for (isolate subset) settings, emit tidy per-gene table

Usage: rev_b2_ctrop_recompute.py cache | agg
"""
import os, gzip, glob, sys, csv
from collections import defaultdict
import numpy as np
from multiprocessing import Pool

BASE = "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
CAN = os.path.join(BASE, "analysis/candida")
REF = os.path.join(CAN, "ref/GCF_000006335.3_ASM633v3_genomic.fna")
GFF = os.path.join(CAN, "ref/genomic.gff")
TBL = os.path.join(BASE, "doc/Candida/tables/ctropicalis_analysis_table.csv")
SCRATCH = os.environ.get("REV_B2_SCRATCH",
    "/tmp/claude-1000/-mnt-HuangLab-26T-HuangLab-Data-Tar-Nine-Ta/459df96a-3261-4f14-8f94-382b2767e924/scratchpad")
CACHE = os.path.join(SCRATCH, "rev_b2_ctrop_cache.npz")
NPROC = int(os.environ.get("REV_B2_NPROC", "48"))
MIN_QUAL, MIN_DP = 30.0, 10
N_DEPTH_DRAW = int(os.environ.get("REV_B2_CT_DRAWS", "20"))
SEED = 20261012

CODON = {}
_b, _aa, _i = "TCAG", "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", 0
for a in _b:
    for b in _b:
        for c in _b:
            CODON[a + b + c] = _aa[_i]; _i += 1
COMP = str.maketrans("ACGTN", "TGCAN")
assert CODON["CTG"] == "S" and CODON["CTA"] == "L" and CODON["ATG"] == "M"
ALT_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
ALT_CHR = "ACGT"


def load_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


site_cache = {}


def syn_sites(codon):
    if codon in site_cache:
        return site_cache[codon]
    aa = CODON.get(codon)
    if aa is None or aa == "*":
        site_cache[codon] = 0.0
        return 0.0
    s = 0.0
    for p in range(3):
        for b in "ACGT":
            if b == codon[p]:
                continue
            if CODON.get(codon[:p] + b + codon[p + 1:]) == aa:
                s += 1 / 3
    site_cache[codon] = s
    return s


def isolates():
    iso = []
    with open(TBL) as fh:
        for r in csv.DictReader(fh):
            if str(r.get("qc_pass", "")).upper() == "TRUE":
                iso.append(r["accession"])
    return [a for a in iso if os.path.exists(os.path.join(CAN, "vcf", a + ".vcf.gz"))]


def build_annotation():
    ref = {k: v.upper() for k, v in load_fasta(REF).items()}
    ctg_len = {k: len(v) for k, v in ref.items()}
    genes = {}
    for line in open(GFF):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        pid = at.get("protein_id")
        if not pid:
            continue
        g = genes.setdefault(pid, {"locus": at.get("locus_tag", ""), "ctg": f[0],
                                   "strand": f[6], "exons": [], "product": at.get("product", "")})
        g["exons"].append((int(f[3]) - 1, int(f[4])))
    gene_ids = sorted(genes)
    gidx = {p: i for i, p in enumerate(gene_ids)}
    gene_arr = {c: np.full(l, -1, dtype=np.int32) for c, l in ctg_len.items()}
    cpos_arr = {c: np.full(l, -1, dtype=np.int32) for c, l in ctg_len.items()}
    cds_seq = {}
    for pid in gene_ids:
        g = genes[pid]
        if g["ctg"] not in ref:
            continue
        ex = sorted(g["exons"], key=lambda t: t[0], reverse=(g["strand"] == "-"))
        seq, k, c = [], 0, g["ctg"]
        for s, e in ex:
            block = ref[c][s:e]
            idxs = np.arange(e - 1, s - 1, -1) if g["strand"] == "-" else np.arange(s, e)
            if g["strand"] == "-":
                block = block.translate(COMP)[::-1]
            seq.append(block)
            gene_arr[c][idxs] = gidx[pid]
            cpos_arr[c][idxs] = np.arange(k, k + len(idxs))
            k += len(idxs)
        cds_seq[pid] = "".join(seq)
    return ctg_len, genes, gene_ids, gidx, gene_arr, cpos_arr, cds_seq


_G = {}


def _vcf_worker(acc):
    gene_arr = _G["gene_arr"]; ctg_i = _G["ctg_i"]
    keys, dose = [], []
    with gzip.open(os.path.join(CAN, "vcf", acc + ".vcf.gz"), "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            r, a = f[3], f[4]
            if len(r) != 1 or "," in a or len(a) != 1:
                continue
            if r not in ALT_CODE or a not in ALT_CODE:
                continue
            try:
                if float(f[5]) < MIN_QUAL:
                    continue
            except ValueError:
                continue
            fmt = f[8].split(":"); smp = f[9].split(":")
            gt = smp[0].replace("|", "/")
            if gt in ("0/0", "./."):
                continue
            if "DP" in fmt:
                try:
                    if int(smp[fmt.index("DP")]) < MIN_DP:
                        continue
                except (ValueError, IndexError):
                    pass
            ci = ctg_i.get(f[0])
            if ci is None:
                continue
            p0 = int(f[1]) - 1
            ga = gene_arr[f[0]]
            if p0 >= ga.shape[0] or ga[p0] < 0:
                continue
            keys.append((ci << 26) | (p0 << 2) | ALT_CODE[a])
            dose.append(2 if gt == "1/1" else 1)
    k = np.asarray(keys, dtype=np.int64)
    d = np.asarray(dose, dtype=np.int8)
    if k.size:
        o = np.argsort(k, kind="stable")
        k, d = k[o], d[o]
        uk, first = np.unique(k, return_index=True)
        d = d[first]
        k = uk
    return k, d


def do_cache():
    print("annotation...", flush=True)
    ctg_len, genes, gene_ids, gidx, gene_arr, cpos_arr, cds_seq = build_annotation()
    ngene = len(gene_ids)
    ctg_order = sorted(ctg_len)
    print(f"  {ngene} transcripts", flush=True)

    Ssit = np.zeros(ngene); Nsit = np.zeros(ngene)
    for pid, s in cds_seq.items():
        S, n = 0.0, 0
        for i in range(0, len(s) - 2, 3):
            cod = s[i:i + 3]
            if len(cod) < 3 or "N" in cod:
                continue
            S += syn_sites(cod); n += 1
        gi = gidx[pid]
        Ssit[gi] = S; Nsit[gi] = 3 * n - S

    iso = isolates()
    print(f"{len(iso)} QC-pass isolates with VCF", flush=True)
    _G["gene_arr"] = gene_arr
    _G["ctg_i"] = {c: i for i, c in enumerate(ctg_order)}
    with Pool(NPROC) as pool:
        res = pool.map(_vcf_worker, iso, chunksize=1)

    allkeys = np.unique(np.concatenate([k for k, d in res]))
    nv = allkeys.size
    print(f"  {nv} unique CDS SNV alleles", flush=True)
    indptr = np.zeros(len(iso) + 1, dtype=np.int64)
    ip, dp = [], []
    for i, (k, d) in enumerate(res):
        vi = np.searchsorted(allkeys, k)
        ip.append(vi.astype(np.int32)); dp.append(d)
        indptr[i + 1] = indptr[i] + vi.size
    vidx = np.concatenate(ip).astype(np.int32)
    vdose = np.concatenate(dp).astype(np.int8)

    print("classifying...", flush=True)
    v_gene = np.full(nv, -1, dtype=np.int32)
    v_class = np.zeros(nv, dtype=np.int8)
    for vi in range(nv):
        key = int(allkeys[vi])
        ci = key >> 26
        p0 = (key >> 2) & ((1 << 24) - 1)
        alt = ALT_CHR[key & 3]
        c = ctg_order[ci]
        gi = int(gene_arr[c][p0])
        if gi < 0:
            continue
        pid = gene_ids[gi]
        s = cds_seq.get(pid)
        if s is None:
            continue
        cpi = int(cpos_arr[c][p0]); codon_i, cp = cpi // 3, cpi % 3
        st = codon_i * 3
        if st + 3 > len(s):
            continue
        cod = s[st:st + 3]
        if "N" in cod:
            continue
        a = alt if genes[pid]["strand"] == "+" else alt.translate(COMP)
        new = cod[:cp] + a + cod[cp + 1:]
        aa0, aa1 = CODON.get(cod), CODON.get(new)
        if aa0 is None or aa1 is None or aa0 == "*":
            continue
        v_gene[vi] = gi
        v_class[vi] = 0 if aa0 == aa1 else (2 if aa1 == "*" else 1)

    # callability from mosdepth per-gene regions (locus_tag keyed)
    print("callability...", flush=True)
    depth_ok = defaultdict(int)
    for acc in iso:
        p = os.path.join(CAN, "depth", acc + ".regions.bed.gz")
        if not os.path.exists(p):
            continue
        with gzip.open(p, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                try:
                    if float(f[4]) >= 10:
                        depth_ok[f[3]] += 1
                except ValueError:
                    continue
    ncall = np.array([depth_ok.get(genes[p]["locus"], 0) for p in gene_ids], dtype=np.int32)

    np.savez_compressed(CACHE, iso=np.array(iso), gene_ids=np.array(gene_ids),
                        locus=np.array([genes[p]["locus"] for p in gene_ids]),
                        product=np.array([genes[p]["product"] for p in gene_ids]),
                        prot_len=np.array([len(cds_seq[p]) // 3 for p in gene_ids], dtype=np.int32),
                        syn_sites=Ssit, nonsyn_sites=Nsit,
                        vidx=vidx, vdose=vdose, indptr=indptr,
                        v_gene=v_gene, v_class=v_class, n_iso_callable=ncall)
    print("wrote", CACHE, flush=True)


def aggregate(z, sample_idx, mac=1):
    vidx, vdose, indptr = z["vidx"], z["vdose"], z["indptr"]
    v_gene, v_class = z["v_gene"], z["v_class"]
    ngene = z["gene_ids"].size
    nv = v_gene.size
    n = 2 * len(sample_idx)
    cnt = np.zeros(nv, dtype=np.int32)
    for g in sample_idx:
        a, b = indptr[g], indptr[g + 1]
        np.add.at(cnt, vidx[a:b], vdose[a:b].astype(np.int32))
    keep = (v_gene >= 0) & (cnt >= mac)
    k = cnt[keep].astype(np.float64)
    gi = v_gene[keep]; cl = v_class[keep]
    pi = 2.0 * k * (n - k) / (n * (n - 1)) if n > 1 else np.zeros_like(k)
    syn = cl == 0; non = ~syn; stop = cl == 2
    return dict(
        syn_seg=np.bincount(gi[syn], minlength=ngene).astype(np.int32),
        nonsyn_seg=np.bincount(gi[non], minlength=ngene).astype(np.int32),
        stop_seg=np.bincount(gi[stop], minlength=ngene).astype(np.int32),
        pi_syn_sum=np.bincount(gi[syn], weights=pi[syn], minlength=ngene),
        pi_nonsyn_sum=np.bincount(gi[non], weights=pi[non], minlength=ngene),
        n_iso=len(sample_idx))


def do_agg():
    z = dict(np.load(CACHE, allow_pickle=True))
    iso = list(z["iso"]); aidx = {a: i for i, a in enumerate(iso)}
    ngene = z["gene_ids"].size
    Ssit, Nsit = z["syn_sites"], z["nonsyn_sites"]
    reps = [l.strip() for l in open(os.path.join(CAN, "constraint/ctrop_lineage_representatives.txt")) if l.strip()]
    rep_idx = np.array(sorted(aidx[a] for a in reps if a in aidx))
    all_idx = np.arange(len(iso))
    print(f"{len(iso)} isolates, {len(rep_idx)} lineage representatives", flush=True)

    settings = [("ct_all_mac1", all_idx, 1), ("ct_dedup_mac1", rep_idx, 1),
                ("ct_all_mac2", all_idx, 2), ("ct_dedup_mac2", rep_idx, 2)]
    rng = np.random.default_rng(SEED)
    # sampling-depth matching: shrink the C. tropicalis isolate set so that the
    # median segregating sites per gene approaches the A. fumigatus value
    for nsub in (10, 20, 40, 80, 151):
        for r in range(N_DEPTH_DRAW if nsub <= 40 else 10):
            pick = rng.choice(rep_idx, size=min(nsub, len(rep_idx)), replace=False)
            settings.append((f"ct_dedupsub{nsub:03d}_r{r:02d}", np.sort(pick), 1))

    gfile = os.path.join(CAN, "constraint/rev_b2_ctrop_gene_static.tsv")
    with open(gfile, "w") as fh:
        fh.write("protein_id\tlocus_tag\tprotein_len\tsyn_sites\tnonsyn_sites\tn_iso_callable\tproduct\n")
        for i in range(ngene):
            fh.write(f"{z['gene_ids'][i]}\t{z['locus'][i]}\t{z['prot_len'][i]}\t{Ssit[i]:.1f}\t"
                     f"{Nsit[i]:.1f}\t{z['n_iso_callable'][i]}\t{z['product'][i]}\n")
    outp = os.path.join(CAN, "constraint/rev_b2_ctrop_setting_counts.tsv.gz")
    with gzip.open(outp, "wt") as fh:
        fh.write("setting\tn_iso\tmac\tprotein_id\tsyn_seg\tnonsyn_seg\tstop_seg\tpi_syn\tpi_nonsyn\n")
        for j, (nm, idx, mac) in enumerate(settings):
            a = aggregate(z, idx, mac)
            piS = np.divide(a["pi_syn_sum"], Ssit, out=np.full(ngene, np.nan), where=Ssit > 0)
            piN = np.divide(a["pi_nonsyn_sum"], Nsit, out=np.full(ngene, np.nan), where=Nsit > 0)
            for i in range(ngene):
                fh.write(f"{nm}\t{a['n_iso']}\t{mac}\t{z['gene_ids'][i]}\t{a['syn_seg'][i]}\t"
                         f"{a['nonsyn_seg'][i]}\t{a['stop_seg'][i]}\t{piS[i]:.6g}\t{piN[i]:.6g}\n")
            if (j + 1) % 10 == 0 or j == 0:
                print(f"  {j+1}/{len(settings)} {nm}", flush=True)
    print("wrote", outp, flush=True)


def do_agg_split():
    """Disjoint split-half replicates of the lineage-representative set, for
    within-genus estimator reliability (see the Aspergillus counterpart)."""
    z = dict(np.load(CACHE, allow_pickle=True))
    iso = list(z["iso"]); aidx = {a: i for i, a in enumerate(iso)}
    ngene = z["gene_ids"].size
    Ssit, Nsit = z["syn_sites"], z["nonsyn_sites"]
    reps = [l.strip() for l in open(os.path.join(CAN, "constraint/ctrop_lineage_representatives.txt")) if l.strip()]
    rep_idx = np.array(sorted(aidx[a] for a in reps if a in aidx))
    rng = np.random.default_rng(SEED + 7)
    ndraw = int(os.environ.get("REV_B2_NSPLIT", "20"))
    outp = os.path.join(CAN, "constraint/rev_b2_ctrop_splithalf_counts.tsv.gz")
    with gzip.open(outp, "wt") as fh:
        fh.write("setting\tn_iso\tmac\tprotein_id\tsyn_seg\tnonsyn_seg\tstop_seg\tpi_syn\tpi_nonsyn\n")
        for r in range(ndraw):
            perm = rng.permutation(rep_idx)
            for h, idx in (("A", np.sort(perm[:len(perm) // 2])), ("B", np.sort(perm[len(perm) // 2:]))):
                a = aggregate(z, idx, 1)
                piS = np.divide(a["pi_syn_sum"], Ssit, out=np.full(ngene, np.nan), where=Ssit > 0)
                piN = np.divide(a["pi_nonsyn_sum"], Nsit, out=np.full(ngene, np.nan), where=Nsit > 0)
                nm = f"ctsplit{r:02d}{h}"
                for i in range(ngene):
                    fh.write(f"{nm}\t{a['n_iso']}\t1\t{z['gene_ids'][i]}\t{a['syn_seg'][i]}\t"
                             f"{a['nonsyn_seg'][i]}\t{a['stop_seg'][i]}\t{piS[i]:.6g}\t{piN[i]:.6g}\n")
            print(f"  split {r+1}/{ndraw}", flush=True)
    print("wrote", outp, flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "cache"
    if cmd == "cache":
        do_cache()
    elif cmd == "agg":
        do_agg()
    elif cmd == "agg_split":
        do_agg_split()
    else:
        sys.exit("unknown " + cmd)
