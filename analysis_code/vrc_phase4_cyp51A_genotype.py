#!/usr/bin/env python3
"""Genotype cyp51A (coding substitutions + promoter tandem repeats TR34/TR46) in each genome.
Non-R step: sequence extraction/alignment/translation has no adequate R equivalent (per lab policy note).
Method: map the Af293 cyp51A locus (gene + ~2 kb promoter) to each assembly with minimap2, extract the
homologous region, align with MAFFT, lift the Af293 CDS through the alignment, translate, and compare.
Promoter insertions are measured in the 500 bp immediately upstream of the start codon."""
import subprocess, sys, os, re
from Bio import SeqIO
from Bio.Seq import Seq

W   = "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta/analysis/vrc_gwas"
C   = f"{W}/04_association/cyp51A"
LOC = f"{C}/cyp51A_locus_Af293.fa"
# Af293 locus was cut as NC_007197.1:1777000-1783900 (1-based). cyp51A is on the minus strand;
# CDS = complement(join(1780204..1781550, 1781622..1781822)).
OFF = 1777000
EXONS = [(1780204,1781550),(1781622,1781822)]     # genome coords, ascending
def loc_coords(a,b): return (a-OFF, b-OFF)        # 0-based within locus

def run(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)

ref = str(next(SeqIO.parse(LOC,"fasta")).seq).upper()

def extract_region(genome, out):
    """minimap2 the locus onto the assembly, take the best hit, extract with flanks."""
    p = run(f"minimap2 -cx asm5 --secondary=no {genome} {LOC} 2>/dev/null")
    best = None
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12: continue
        block = int(f[10])
        if best is None or block > best[0]: best = (block, f)
    if best is None: return None
    f = best[1]; tname, tstart, tend, strand = f[5], int(f[7]), int(f[8]), f[4]
    tlen = int(f[6])
    s = max(0, tstart-1500); e = min(tlen, tend+1500)
    r = run(f"samtools faidx {genome} {tname}:{s+1}-{e}")
    if r.returncode: return None
    seq = "".join(l.strip() for l in r.stdout.splitlines() if not l.startswith(">"))
    if strand == "-": seq = str(Seq(seq).reverse_complement())
    open(out,"w").write(f">query\n{seq}\n")
    return dict(contig=tname, start=s+1, end=e, strand=strand, matchlen=best[0])

def align(qfile, out):
    run(f"cat {LOC} {qfile} > {out}.in.fa")
    r = run(f"mafft --quiet --auto --thread 4 {out}.in.fa > {out}")
    return r.returncode == 0

def lift(aln_ref, aln_qry, pos0):
    """map a 0-based ungapped ref position to the aligned-column index"""
    n = -1
    for i,ch in enumerate(aln_ref):
        if ch != "-":
            n += 1
            if n == pos0: return i
    return None

def genotype(acc, genome):
    q = f"{C}/tmp_{acc}.fa"; a = f"{C}/tmp_{acc}.aln"
    info = extract_region(genome, q)
    if info is None: return dict(accession=acc, status="no cyp51A hit")
    if not align(q, a): return dict(accession=acc, status="alignment failed")
    recs = list(SeqIO.parse(a,"fasta"))
    ar = str([r for r in recs if r.id!="query"][0].seq).upper()
    aq = str([r for r in recs if r.id=="query"][0].seq).upper()
    # --- CDS: build ref & query CDS by lifting exon coordinates
    cds_r, cds_q = [], []
    for (g1,g2) in EXONS:
        s0,e0 = loc_coords(g1,g2)
        ci, cj = lift(ar,aq,s0), lift(ar,aq,e0)
        if ci is None or cj is None: return dict(accession=acc, status="CDS lift failed")
        cds_r.append(ar[ci:cj+1]); cds_q.append(aq[ci:cj+1])
    # exons are ascending in genome coords but gene is minus strand -> revcomp the concatenation
    rr = Seq("".join(cds_r).replace("-","")).reverse_complement()
    qq_gapped = "".join(cds_q)
    qq = Seq(qq_gapped.replace("-","")).reverse_complement()
    prot_r = str(rr.translate()); prot_q = str(qq.translate())
    subs = []
    if len(prot_q) == len(prot_r):
        subs = [f"{prot_r[i]}{i+1}{prot_q[i]}" for i in range(len(prot_r)) if prot_r[i]!=prot_q[i]]
        frameshift = ""
    else:
        frameshift = f"length {len(prot_q)} vs {len(prot_r)}"
    # --- promoter: 500 bp immediately upstream of ATG (minus strand -> higher genome coords)
    atg0 = loc_coords(EXONS[1][1], EXONS[1][1])[0]          # position of the A of ATG in locus coords
    pi, pj = lift(ar,aq,atg0), lift(ar,aq,min(atg0+500, len(ref)-1))
    ins = 0; prom_r = prom_q = ""
    if pi is not None and pj is not None:
        prom_r = ar[pi:pj+1]; prom_q = aq[pi:pj+1]
        ins = len(prom_q.replace("-","")) - len(prom_r.replace("-",""))
    tr = "TR34" if 30 <= ins <= 38 else "TR46" if 42 <= ins <= 50 else \
         "TR34/3" if 100 <= ins <= 108 else "TR34/2" if 66 <= ins <= 74 else \
         ("none" if abs(ins) <= 3 else f"other_ins_{ins}bp")
    for f in (q, a, a+".in.fa"):
        try: os.remove(f)
        except OSError: pass
    return dict(accession=acc, status="ok", contig=info["contig"], promoter_insert_bp=ins,
                promoter_genotype=tr, protein_len=len(prot_q), frameshift=frameshift,
                aa_substitutions=";".join(subs), n_aa_subs=len(subs))

if __name__ == "__main__":
    import csv, glob
    genomes = sorted(glob.glob(f"{W}/genomes/*.fna"))
    rows = []
    for i,g in enumerate(genomes,1):
        acc = os.path.basename(g)[:-4]
        try: r = genotype(acc, os.path.realpath(g))
        except Exception as e: r = dict(accession=acc, status=f"error: {e}")
        rows.append(r); print(f"[{i}/{len(genomes)}] {acc} {r.get('promoter_genotype','-')} {r.get('aa_substitutions','')[:60]}", flush=True)
    cols = ["accession","status","contig","promoter_insert_bp","promoter_genotype","protein_len","frameshift","n_aa_subs","aa_substitutions"]
    with open(f"{C}/cyp51A_genotypes.csv","w",newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    print("wrote", f"{C}/cyp51A_genotypes.csv")
