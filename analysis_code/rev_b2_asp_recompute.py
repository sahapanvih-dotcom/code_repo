#!/usr/bin/env python3
"""rev_b2 (manuscript revision, analysis B2): recompute A. fumigatus per-gene piN/piS
under many analytical settings and under repeated balanced lineage subsampling.

Why Python and not R: this reuses the project's existing variant parser / Nei-Gojobori
site-counting / codon-classification implementation from analysis_code/td_constraint.py
verbatim (same GFF handling, same CDS build, same syn_sites(), same
pi = 2k(n-k)/(n(n-1)) estimator, same callability definition). Re-implementing that
parser in R would be wasteful and would risk introducing differences that would be
mistaken for instability. ALL summarising, statistics and plotting is done in R
(analysis_code/rev_b2_analyse.R, rev_b2_crossgenus.R, rev_b2_figure.R).

Two stages:
  cache : parse all 315 genomes once -> compact per-genome variant index + per-gene
          covered-bp matrix, cached as .npz
  agg   : aggregate the cache for an arbitrary list of (genome subset x allele-count
          threshold) settings and emit tidy per-gene count tables

Usage:
  rev_b2_asp_recompute.py cache
  rev_b2_asp_recompute.py agg
"""
import os, gzip, glob, sys, json
from collections import defaultdict
import numpy as np
from multiprocessing import Pool

BASE = "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
TD = os.path.join(BASE, "analysis/target_discovery")
REF = os.path.join(TD, "ref/Af293.fna")
GFF = os.path.join(TD, "ref/Af293.gff")
SCRATCH = os.environ.get("REV_B2_SCRATCH",
    "/tmp/claude-1000/-mnt-HuangLab-26T-HuangLab-Data-Tar-Nine-Ta/459df96a-3261-4f14-8f94-382b2767e924/scratchpad")
os.makedirs(SCRATCH, exist_ok=True)
CACHE = os.path.join(SCRATCH, "rev_b2_asp_cache.npz")
NPROC = int(os.environ.get("REV_B2_NPROC", "48"))

N_RESAMPLE = int(os.environ.get("REV_B2_NRESAMPLE", "100"))
RESAMPLE_SEED = 20261012

# ---------------------------------------------------------------- genetic code
CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in _b:
    for b in _b:
        for c in _b:
            CODON[a + b + c] = _aa[i]; i += 1
COMP = str.maketrans("ACGTN", "TGCAN")
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


# ---------------------------------------------------------------- annotation
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
                                   "strand": f[6], "exons": [],
                                   "product": at.get("product", "")})
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
    return ref, ctg_len, genes, gene_ids, gidx, gene_arr, cpos_arr, cds_seq


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


# ---------------------------------------------------------------- workers
_G = {}


def _init_worker(gene_arr, ctg_order, ctg_len):
    _G["gene_arr"] = gene_arr
    _G["ctg_i"] = {c: i for i, c in enumerate(ctg_order)}
    _G["ctg_len"] = ctg_len


def _snp_worker(path):
    """Return (encoded_variant_keys, per_genome_multiplicity) for CDS variants."""
    gene_arr = _G["gene_arr"]; ctg_i = _G["ctg_i"]
    keys = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            c, p, r, a = f[0], f[1], f[2], f[3]
            if a not in ALT_CODE or r not in ALT_CODE:
                continue
            ci = ctg_i.get(c)
            if ci is None:
                continue
            pos = int(p)
            ga = gene_arr[c]
            if pos >= ga.shape[0] or ga[pos] < 0:
                continue
            keys.append((ci << 26) | (pos << 2) | ALT_CODE[a])
    arr = np.asarray(keys, dtype=np.int64)
    u, m = np.unique(arr, return_counts=True)
    return u, m.astype(np.int16)


def _cov_worker(path):
    gene_arr = _G["gene_arr"]; ctg_len = _G["ctg_len"]
    ngene = _G["ngene"]
    out = np.zeros(ngene, dtype=np.int32)
    seen = {}
    for line in open(path):
        f = line.split()
        if len(f) < 3:
            continue
        c = f[0]
        if c not in ctg_len:
            continue
        if c not in seen:
            seen[c] = np.zeros(ctg_len[c], dtype=bool)
        seen[c][int(f[1]):int(f[2])] = True
    for c, mask in seen.items():
        ga = gene_arr[c]
        sel = ga[mask]
        sel = sel[sel >= 0]
        if sel.size:
            out += np.bincount(sel, minlength=ngene).astype(np.int32)
    return out


def _init_worker2(gene_arr, ctg_order, ctg_len, ngene):
    _init_worker(gene_arr, ctg_order, ctg_len)
    _G["ngene"] = ngene


# ---------------------------------------------------------------- cache stage
def do_cache():
    print("building annotation...", flush=True)
    ref, ctg_len, genes, gene_ids, gidx, gene_arr, cpos_arr, cds_seq = build_annotation()
    ngene = len(gene_ids)
    print(f"  {ngene} transcripts", flush=True)

    # Nei-Gojobori site counts
    print("counting syn/nonsyn sites...", flush=True)
    Ssit = np.zeros(ngene); Nsit = np.zeros(ngene); ncod = np.zeros(ngene, dtype=np.int32)
    for pid, s in cds_seq.items():
        S = 0.0; n = 0
        for i in range(0, len(s) - 2, 3):
            cod = s[i:i + 3]
            if len(cod) < 3 or "N" in cod:
                continue
            S += syn_sites(cod); n += 1
        gi = gidx[pid]
        Ssit[gi] = S; Nsit[gi] = 3 * n - S; ncod[gi] = n

    # CDS bp per gene actually owned in gene_arr (matches td_constraint.py's coverage denominator)
    own_bp = np.zeros(ngene, dtype=np.int32)
    for c, ga in gene_arr.items():
        sel = ga[ga >= 0]
        if sel.size:
            own_bp += np.bincount(sel, minlength=ngene).astype(np.int32)

    ctg_order = sorted(ctg_len)
    snp_files = sorted(glob.glob(os.path.join(TD, "var/*.snp.tsv.gz")))
    accs = [os.path.basename(f).replace(".snp.tsv.gz", "") for f in snp_files]
    cov_files = [os.path.join(TD, "var", a + ".cov.bed") for a in accs]
    assert all(os.path.exists(f) for f in cov_files), "missing cov.bed"
    ngen = len(accs)
    print(f"{ngen} genomes", flush=True)

    # populate the worker globals in the PARENT so forked workers inherit them
    # copy-on-write instead of pickling ~230 MB of coordinate arrays per worker
    _init_worker2(gene_arr, ctg_order, ctg_len, ngene)

    print("parsing SNP files...", flush=True)
    with Pool(NPROC) as pool:
        snp_res = pool.map(_snp_worker, snp_files, chunksize=1)
        print("parsing coverage beds...", flush=True)
        cov_res = pool.map(_cov_worker, cov_files, chunksize=1)

    allkeys = np.unique(np.concatenate([u for u, m in snp_res]))
    nv = allkeys.size
    print(f"  {nv} unique CDS variant alleles", flush=True)

    # per-genome sparse index into allkeys
    indptr = np.zeros(ngen + 1, dtype=np.int64)
    idx_parts, mult_parts = [], []
    for i, (u, m) in enumerate(snp_res):
        vi = np.searchsorted(allkeys, u)
        idx_parts.append(vi.astype(np.int32)); mult_parts.append(m)
        indptr[i + 1] = indptr[i] + vi.size
    vidx = np.concatenate(idx_parts).astype(np.int32)
    vmult = np.concatenate(mult_parts).astype(np.int16)
    print(f"  {vidx.size} genome-variant observations "
          f"({int((vmult > 1).sum())} duplicated alignment-block lines)", flush=True)

    covmat = np.vstack(cov_res)  # ngen x ngene covered CDS bp

    # ---- classify every unique variant once -------------------------------
    print("classifying variants...", flush=True)
    v_gene = np.full(nv, -1, dtype=np.int32)
    v_class = np.zeros(nv, dtype=np.int8)   # 0 syn, 1 nonsyn (non-stop), 2 stop-gained
    for vi in range(nv):
        key = int(allkeys[vi])
        ci = key >> 26
        pos = (key >> 2) & ((1 << 24) - 1)
        alt = ALT_CHR[key & 3]
        c = ctg_order[ci]
        gi = int(gene_arr[c][pos])
        if gi < 0:
            continue
        pid = gene_ids[gi]
        s = cds_seq.get(pid)
        if s is None:
            continue
        cpi = int(cpos_arr[c][pos])
        codon_i, cp = cpi // 3, cpi % 3
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
        if aa0 == aa1:
            v_class[vi] = 0
        elif aa1 == "*":
            v_class[vi] = 2
        else:
            v_class[vi] = 1
    print(f"  {int((v_gene >= 0).sum())} classifiable", flush=True)

    np.savez_compressed(
        CACHE,
        accs=np.array(accs), gene_ids=np.array(gene_ids),
        locus=np.array([genes[p]["locus"] for p in gene_ids]),
        product=np.array([genes[p]["product"] for p in gene_ids]),
        prot_len=np.array([len(cds_seq[p]) // 3 if p in cds_seq else 0 for p in gene_ids], dtype=np.int32),
        n_codons=ncod, syn_sites=Ssit, nonsyn_sites=Nsit, own_bp=own_bp,
        vidx=vidx, vmult=vmult, indptr=indptr,
        v_gene=v_gene, v_class=v_class, covmat=covmat)
    print("wrote", CACHE, flush=True)


# ---------------------------------------------------------------- agg stage
def aggregate(z, sample_idx, mac, dedup_lines=True):
    """Per-gene counts/pi for a genome subset at derived-allele-count threshold `mac`."""
    vidx, vmult, indptr = z["vidx"], z["vmult"], z["indptr"]
    v_gene, v_class = z["v_gene"], z["v_class"]
    ngene = z["gene_ids"].size
    nv = v_gene.size
    n = len(sample_idx)
    cnt = np.zeros(nv, dtype=np.int32)
    for g in sample_idx:
        a, b = indptr[g], indptr[g + 1]
        sl = vidx[a:b]
        if dedup_lines:
            np.add.at(cnt, sl, 1)
        else:
            np.add.at(cnt, sl, vmult[a:b].astype(np.int32))
    keep = (v_gene >= 0) & (cnt >= mac)
    k = cnt[keep].astype(np.float64)
    gi = v_gene[keep]
    cl = v_class[keep]
    pi = 2.0 * k * (n - k) / (n * (n - 1)) if n > 1 else np.zeros_like(k)
    syn = cl == 0
    non = ~syn
    stop = cl == 2
    out = {}
    out["syn_seg"] = np.bincount(gi[syn], minlength=ngene).astype(np.int32)
    out["nonsyn_seg"] = np.bincount(gi[non], minlength=ngene).astype(np.int32)
    out["stop_seg"] = np.bincount(gi[stop], minlength=ngene).astype(np.int32)
    out["pi_syn_sum"] = np.bincount(gi[syn], weights=pi[syn], minlength=ngene)
    out["pi_nonsyn_sum"] = np.bincount(gi[non], weights=pi[non], minlength=ngene)
    cov = z["covmat"][sample_idx, :].sum(axis=0).astype(np.float64)
    own = z["own_bp"].astype(np.float64)
    out["mean_cov_genomes"] = np.divide(cov, own, out=np.zeros(ngene), where=own > 0)
    out["n_genomes"] = n
    return out


def do_agg():
    z = dict(np.load(CACHE, allow_pickle=True))
    accs = list(z["accs"])
    aidx = {a: i for i, a in enumerate(accs)}
    gene_ids = z["gene_ids"]; ngene = gene_ids.size
    Ssit, Nsit = z["syn_sites"], z["nonsyn_sites"]

    reps = [l.strip() for l in open(os.path.join(TD, "popstruct/representatives.txt")) if l.strip()]
    rep_idx = np.array(sorted(aidx[a] for a in reps if a in aidx))
    all_idx = np.arange(len(accs))

    # lineage membership
    lin = defaultdict(list)
    with open(os.path.join(TD, "popstruct/lineages.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in aidx:
                lin[int(f[1])].append(aidx[f[0]])
    lineages = [np.array(v) for k, v in sorted(lin.items())]
    print(f"{len(accs)} genomes, {len(rep_idx)} representatives, {len(lineages)} lineages", flush=True)

    settings = []
    # --- published and exact-reproduction checks
    settings.append(("published_repro_perline", all_idx, 1, False))
    settings.append(("all315_mac1", all_idx, 1, True))
    # --- one-factor: derived-allele-count / frequency threshold, all 315
    for mac, nm in [(2, "all315_mac2"), (3, "all315_mac3"), (5, "all315_mac5"),
                    (4, "all315_maf0.01"), (16, "all315_maf0.05")]:
        settings.append((nm, all_idx, mac, True))
    # --- one-factor: lineage de-duplication (published secondary analysis)
    settings.append(("dedup188_mac1", rep_idx, 1, True))
    # --- combination: de-duplication x allele-count threshold
    for mac, nm in [(2, "dedup188_mac2"), (3, "dedup188_mac3"), (10, "dedup188_maf0.05")]:
        settings.append((nm, rep_idx, mac, True))

    rng = np.random.default_rng(RESAMPLE_SEED)
    for r in range(N_RESAMPLE):
        pick = np.array([lg[rng.integers(len(lg))] for lg in lineages])
        settings.append((f"resample{r:03d}", pick, 1, True))

    outp = os.path.join(BASE, "analysis/target_discovery/rev_b2_setting_counts.tsv.gz")
    gfile = os.path.join(BASE, "analysis/target_discovery/rev_b2_gene_static.tsv")
    with open(gfile, "w") as fh:
        fh.write("protein_id\tlocus_tag\tprotein_len\tn_codons\tsyn_sites\tnonsyn_sites\tproduct\n")
        for i in range(ngene):
            fh.write(f"{gene_ids[i]}\t{z['locus'][i]}\t{z['prot_len'][i]}\t{z['n_codons'][i]}\t"
                     f"{Ssit[i]:.1f}\t{Nsit[i]:.1f}\t{z['product'][i]}\n")
    print("wrote", gfile, flush=True)

    with gzip.open(outp, "wt") as fh:
        fh.write("setting\tn_genomes\tmac\tprotein_id\tsyn_seg\tnonsyn_seg\tstop_seg\t"
                 "pi_syn\tpi_nonsyn\tmean_cov_genomes\n")
        for j, (nm, idx, mac, dd) in enumerate(settings):
            a = aggregate(z, idx, mac, dedup_lines=dd)
            piS = np.divide(a["pi_syn_sum"], Ssit, out=np.full(ngene, np.nan), where=Ssit > 0)
            piN = np.divide(a["pi_nonsyn_sum"], Nsit, out=np.full(ngene, np.nan), where=Nsit > 0)
            n = a["n_genomes"]
            for i in range(ngene):
                if Ssit[i] <= 0 and Nsit[i] <= 0:
                    continue
                fh.write(f"{nm}\t{n}\t{mac}\t{gene_ids[i]}\t{a['syn_seg'][i]}\t"
                         f"{a['nonsyn_seg'][i]}\t{a['stop_seg'][i]}\t{piS[i]:.6g}\t{piN[i]:.6g}\t"
                         f"{a['mean_cov_genomes'][i]:.2f}\n")
            if (j + 1) % 10 == 0 or j == 0:
                print(f"  {j+1}/{len(settings)} {nm} (n={n})", flush=True)
    print("wrote", outp, flush=True)


def do_agg_split():
    """Disjoint split-half replicates, for within-genus reliability of the estimator.

    Two independent halves of the LINEAGE representative set give an upper bound on
    how well any per-gene constraint estimate can correlate with anything: it is the
    repeatability of the estimator itself at this sample size. Needed to interpret
    the cross-genus Spearman rho.
    """
    z = dict(np.load(CACHE, allow_pickle=True))
    accs = list(z["accs"]); aidx = {a: i for i, a in enumerate(accs)}
    gene_ids = z["gene_ids"]; ngene = gene_ids.size
    Ssit, Nsit = z["syn_sites"], z["nonsyn_sites"]
    reps = [l.strip() for l in open(os.path.join(TD, "popstruct/representatives.txt")) if l.strip()]
    rep_idx = np.array(sorted(aidx[a] for a in reps if a in aidx))
    rng = np.random.default_rng(RESAMPLE_SEED + 7)
    ndraw = int(os.environ.get("REV_B2_NSPLIT", "20"))
    outp = os.path.join(BASE, "analysis/target_discovery/rev_b2_splithalf_counts.tsv.gz")
    with gzip.open(outp, "wt") as fh:
        fh.write("setting\tn_genomes\tmac\tprotein_id\tsyn_seg\tnonsyn_seg\tstop_seg\t"
                 "pi_syn\tpi_nonsyn\tmean_cov_genomes\n")
        for r in range(ndraw):
            perm = rng.permutation(rep_idx)
            halves = {"A": np.sort(perm[:len(perm) // 2]), "B": np.sort(perm[len(perm) // 2:])}
            for h, idx in halves.items():
                a = aggregate(z, idx, 1, dedup_lines=True)
                piS = np.divide(a["pi_syn_sum"], Ssit, out=np.full(ngene, np.nan), where=Ssit > 0)
                piN = np.divide(a["pi_nonsyn_sum"], Nsit, out=np.full(ngene, np.nan), where=Nsit > 0)
                nm = f"split{r:02d}{h}"
                for i in range(ngene):
                    fh.write(f"{nm}\t{a['n_genomes']}\t1\t{gene_ids[i]}\t{a['syn_seg'][i]}\t"
                             f"{a['nonsyn_seg'][i]}\t{a['stop_seg'][i]}\t{piS[i]:.6g}\t"
                             f"{piN[i]:.6g}\t{a['mean_cov_genomes'][i]:.2f}\n")
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
        sys.exit("unknown command " + cmd)
