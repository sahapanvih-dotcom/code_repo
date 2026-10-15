#!/usr/bin/env python3
"""Candidate-gene panel: protein-changing variants in known azole-resistance genes, 55 genomes.

Method (the one validated on cyp51A): for each gene take the Af293 genomic locus (CDS span + flank),
map it into the assembly with minimap2, cut the homologous region, align with MAFFT, lift the Af293
exon coordinates through the alignment, translate, and diff against the Af293 protein.
Validation history: this method reproduced 40/41 published cyp51A genotypes, and NCBI's own
annotation of HMR-AF-270 independently confirms its 5-substitution call (DIAMOND 98.1%, 5 mismatches).
miniprot --aln and a cs-tag reconstruction were both tried first and gave wrong answers; discarded."""
import subprocess, os, glob, csv, re, sys
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq

BASE = "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
W = f"{BASE}/analysis/vrc_gwas"; P = f"{W}/04_association/panel"
REF = f"{BASE}/raw_data/coverage_pass/Aspergillus_fumigatus/GCF_000002655.1/GCF_000002655.1_ASM265v1_genomic.fna"
GFF = f"{BASE}/raw_data/coverage_pass/Aspergillus_fumigatus/GCF_000002655.1/genomic.gff"
FAA = f"{BASE}/raw_data/coverage_pass/Aspergillus_fumigatus/GCF_000002655.1/protein.faa"
PANEL_LT = ["AFUA_2G02690","AFUA_2G12910","AFUA_7G04740","AFUA_2G14250","AFUA_3G02340","AFUA_2G10600","AFUA_4G08340"]
FLANK = 500
MMBASE = ("docker run --rm -u {uid}:{gid} -v /mnt/HuangLab_26T:/mnt/HuangLab_26T "
          "quay.io/biocontainers/minimap2:2.31--h118bc1c_0 minimap2 ")

def gene_models():
    """gene -> (chrom, strand, [(s,e) CDS exons ascending], protein_id)"""
    cds = {}
    for line in open(GFF):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        at = dict(kv.split("=",1) for kv in f[8].split(";") if "=" in kv)
        g = at.get("locus_tag")
        if g not in PANEL_LT: continue
        pid = at.get("protein_id","")
        key = (g, pid)
        cds.setdefault(key, dict(chrom=f[0], strand=f[6], exons=[]))["exons"].append((int(f[3]), int(f[4])))
    out = {}
    for (g, pid), d in cds.items():
        if g in out: continue           # first transcript only
        d["exons"].sort(); d["pid"] = pid; out[g] = d
    return out

def run(cmd): return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout

def one_genome(args):
    acc, gpath, models, refprot = args
    rows = []
    for g, mdl in models.items():
        lo = mdl["exons"][0][0] - FLANK; hi = mdl["exons"][-1][1] + FLANK
        loc = f"{P}/tmp2/{acc}_{g}_ref.fa"
        run(f"samtools faidx {REF} {mdl['chrom']}:{max(1,lo)}-{hi} > {loc}")
        paf = run(MMBASE.format(uid=os.getuid(), gid=os.getgid()) + f"-cx asm5 --secondary=no -t 2 {gpath} {loc} 2>/dev/null")
        best = None
        for line in paf.splitlines():
            f = line.split("\t")
            if len(f) < 12: continue
            if best is None or int(f[9]) > best[0]: best = (int(f[9]), f)
        if best is None:
            rows.append(dict(accession=acc, gene=g, status="no_hit", n_subs="", substitutions="")); continue
        f = best[1]; tname, ts, te, strand, tlen = f[5], int(f[7]), int(f[8]), f[4], int(f[6])
        s = max(0, ts-1000); e = min(tlen, te+1000)
        reg = run(f"samtools faidx {gpath} {tname}:{s+1}-{e}")
        seq = "".join(l.strip() for l in reg.splitlines() if not l.startswith(">"))
        if strand == "-": seq = str(Seq(seq).reverse_complement())
        qf = f"{P}/tmp2/{acc}_{g}_q.fa"; open(qf,"w").write(f">query\n{seq}\n")
        af = f"{P}/tmp2/{acc}_{g}.aln"
        run(f"cat {loc} {qf} > {af}.in && mafft --quiet --auto --thread 1 {af}.in > {af}")
        try:
            recs = list(SeqIO.parse(af,"fasta"))
            ar = str([r for r in recs if r.id!="query"][0].seq).upper()
            aq = str([r for r in recs if r.id=="query"][0].seq).upper()
        except Exception:
            rows.append(dict(accession=acc, gene=g, status="aln_fail", n_subs="", substitutions="")); continue
        # map ungapped ref offsets -> alignment columns
        col, n = {}, -1
        for i,ch in enumerate(ar):
            if ch != "-":
                n += 1; col[n] = i
        off = max(1, lo)
        rr, qq = [], []
        ok = True
        for (gs, ge) in mdl["exons"]:
            a, b = gs-off, ge-off
            if a not in col or b not in col: ok = False; break
            rr.append(ar[col[a]:col[b]+1]); qq.append(aq[col[a]:col[b]+1])
        for x in (loc, qf, af, af+".in"):
            try: os.remove(x)
            except OSError: pass
        if not ok:
            rows.append(dict(accession=acc, gene=g, status="lift_fail", n_subs="", substitutions="")); continue
        rs = Seq("".join(rr).replace("-","")); qs_ = Seq("".join(qq).replace("-",""))
        if mdl["strand"] == "-": rs, qs_ = rs.reverse_complement(), qs_.reverse_complement()
        pr = str(rs.translate()).rstrip("*"); pq = str(qs_.translate()).rstrip("*")
        if pr != refprot[g]:
            rows.append(dict(accession=acc, gene=g, status="ref_model_mismatch", n_subs="", substitutions="")); continue
        if len(pq) != len(pr):
            rows.append(dict(accession=acc, gene=g, status=f"length_diff({len(pq)}v{len(pr)})", n_subs="", substitutions="")); continue
        subs = [f"{pr[i]}{i+1}{pq[i]}" for i in range(len(pr)) if pr[i]!=pq[i] and pq[i] not in "X*"]
        rows.append(dict(accession=acc, gene=g, status="ok", n_subs=len(subs), substitutions=";".join(subs)))
    return rows

def main():
    os.makedirs(f"{P}/tmp2", exist_ok=True)
    models = gene_models()
    prot = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(FAA,"fasta")}
    refprot = {g: prot.get(m["pid"],"") for g,m in models.items()}
    models = {g:m for g,m in models.items() if refprot.get(g)}
    print(f"extended panel gene models: {len(models)} -> {sorted(models)}", flush=True)
    genomes = sorted(glob.glob(f"{W}/02_features/genomes_real/*.fna"))
    jobs = [(os.path.basename(p)[:-4], p, models, refprot) for p in genomes]
    allrows = []
    with ProcessPoolExecutor(max_workers=14) as ex:
        for i, rows in enumerate(ex.map(one_genome, jobs), 1):
            allrows += rows; print(f"[{i}/{len(jobs)}] {rows[0]['accession']}", flush=True)
    with open(f"{P}/panel_variants_extended.csv","w",newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession","gene","status","n_subs","substitutions"]); w.writeheader(); w.writerows(allrows)
    print("wrote", f"{P}/panel_variants_extended.csv")

if __name__ == "__main__": main()
