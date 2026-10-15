#!/usr/bin/env python3
"""Hunt for tandem-repeat expansions in the promoters of the 10 panel genes across all 55 genomes.

Rationale: the famous azole-resistance mechanism (TR34/TR46) is a duplication in the cyp51A promoter
that raises target expression. Nobody has looked for the same trick at the other panel genes.
Built-in control: the method MUST recover TR34/TR46 at cyp51A, or it is not working.

Method: take 1.5 kb upstream of each gene's CDS start in Af293 (strand-aware), map it into each
assembly with minimap2, extract the homologous region, align with MAFFT, and measure the net
insertion in the 600 bp immediately upstream of the start codon."""
import subprocess, os, glob, csv, sys
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq

BASE = "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"
W = f"{BASE}/analysis/vrc_gwas"; P = f"{W}/06_promoters"
REF = f"{BASE}/raw_data/coverage_pass/Aspergillus_fumigatus/GCF_000002655.1/GCF_000002655.1_ASM265v1_genomic.fna"
GFF = f"{BASE}/raw_data/coverage_pass/Aspergillus_fumigatus/GCF_000002655.1/genomic.gff"
LOCI = {"cyp51A":"AFUA_4G06890","cyp51B":"AFUA_7G03740","hmg1":"AFUA_2G03700","cdr1B":"AFUA_1G14330",
        "atrR":"AFUA_2G02690","ptaB":"AFUA_2G12910","nctA":"AFUA_2G14250","srbA":"AFUA_2G01260",
        "complexI":"AFUA_2G10600","cox10":"AFUA_4G08340","erg25B":"AFUA_4G04820","rttA":"AFUA_7G04740",
        "hapE":"AFUA_6G05300","nctB":"AFUA_3G02340"}
UP = 1500          # promoter length taken from Af293
WIN = 600          # window immediately upstream of ATG in which insertions are measured
MM = ("docker run --rm -u {uid}:{gid} -v /mnt/HuangLab_26T:/mnt/HuangLab_26T "
      "quay.io/biocontainers/minimap2:2.31--h118bc1c_0 minimap2 -cx asm5 --secondary=no -t 2 {g} {q} 2>/dev/null")

def run(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout

def cds_bounds():
    d = {}
    for line in open(GFF):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        at = dict(kv.split("=",1) for kv in f[8].split(";") if "=" in kv)
        lt = at.get("locus_tag")
        if lt not in LOCI.values(): continue
        r = d.setdefault(lt, dict(chrom=f[0], strand=f[6], lo=int(f[3]), hi=int(f[4])))
        r["lo"] = min(r["lo"], int(f[3])); r["hi"] = max(r["hi"], int(f[4]))
    return d

def one(args):
    acc, gpath, prom = args
    rows = []
    for gene, (chrom, strand, s, e, refseq) in prom.items():
        q = f"{P}/tmp/{acc}_{gene}.fa"; open(q,"w").write(f">prom\n{refseq}\n")
        paf = run(MM.format(uid=os.getuid(), gid=os.getgid(), g=gpath, q=q))
        best = None
        for line in paf.splitlines():
            f = line.split("\t")
            if len(f) < 12: continue
            if best is None or int(f[9]) > best[0]: best = (int(f[9]), f)
        if best is None:
            rows.append(dict(accession=acc, gene=gene, status="no_hit", insert_bp="")); continue
        f = best[1]; tn, ts, te, st, tl = f[5], int(f[7]), int(f[8]), f[4], int(f[6])
        a, b = max(0, ts-800), min(tl, te+800)
        reg = run(f"samtools faidx {gpath} {tn}:{a+1}-{b}")
        seq = "".join(l.strip() for l in reg.splitlines() if not l.startswith(">"))
        if st == "-": seq = str(Seq(seq).reverse_complement())
        qf = f"{P}/tmp/{acc}_{gene}_q.fa"; open(qf,"w").write(f">query\n{seq}\n")
        af = f"{P}/tmp/{acc}_{gene}.aln"
        run(f"cat {q} {qf} > {af}.in && mafft --quiet --auto --thread 1 {af}.in > {af}")
        ins = None
        try:
            recs = list(SeqIO.parse(af,"fasta"))
            ar = str([r for r in recs if r.id!="query"][0].seq).upper()
            aq = str([r for r in recs if r.id=="query"][0].seq).upper()
            # ref promoter is oriented so its 3' end abuts the ATG -> measure the last WIN ref bases
            n = 0; start_col = None
            nref = sum(1 for c in ar if c != "-")
            target = max(0, nref - WIN)
            for i, c in enumerate(ar):
                if c != "-":
                    if n == target: start_col = i; break
                    n += 1
            if start_col is not None:
                rr = ar[start_col:].replace("-",""); qq = aq[start_col:].replace("-","")
                ins = len(qq) - len(rr)
        except Exception:
            pass
        for x in (q, qf, af, af+".in"):
            try: os.remove(x)
            except OSError: pass
        rows.append(dict(accession=acc, gene=gene, status="ok" if ins is not None else "aln_fail",
                         insert_bp=ins if ins is not None else ""))
    return rows

def main():
    os.makedirs(f"{P}/tmp", exist_ok=True)
    cb = cds_bounds()
    prom = {}
    for gene, lt in LOCI.items():
        if lt not in cb: print("missing:", gene); continue
        r = cb[lt]
        if r["strand"] == "+": s, e = max(1, r["lo"]-UP), r["lo"]-1
        else:                  s, e = r["hi"]+1, r["hi"]+UP
        seq = run(f"samtools faidx {REF} {r['chrom']}:{s}-{e}")
        seq = "".join(l.strip() for l in seq.splitlines() if not l.startswith(">"))
        if r["strand"] == "-": seq = str(Seq(seq).reverse_complement())
        prom[gene] = (r["chrom"], r["strand"], s, e, seq)
    print(f"promoters extracted: {len(prom)} ({UP} bp each)", flush=True)
    genomes = sorted(glob.glob(f"{W}/02_features/genomes_real/*.fna"))
    jobs = [(os.path.basename(p)[:-4], p, prom) for p in genomes]
    allrows = []
    with ProcessPoolExecutor(max_workers=14) as ex:
        for i, rows in enumerate(ex.map(one, jobs), 1):
            allrows += rows; print(f"[{i}/{len(jobs)}] {rows[0]['accession']}", flush=True)
    with open(f"{P}/promoter_insertions.csv","w",newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession","gene","status","insert_bp"]); w.writeheader(); w.writerows(allrows)
    print("wrote", f"{P}/promoter_insertions.csv")

if __name__ == "__main__": main()
