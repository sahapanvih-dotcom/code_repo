#!/usr/bin/env python3
"""Call Erg11p substitutions per isolate, with zygosity.
C. tropicalis is diploid, so each isolate's two alleles are reconstructed separately
(bcftools consensus -H 1 / -H 2) and translated. A substitution present on both alleles is
homozygous, on one allele heterozygous - matching how Fan 2023 reports its genotypes."""
import subprocess, os, glob, csv, sys
from Bio import SeqIO
from Bio.Seq import Seq
C = "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta/analysis/candida"
REG = f"{C}/erg11/region.fa"
BCF = ("docker run --rm -u %d:%d -v /mnt/HuangLab_26T:/mnt/HuangLab_26T -v /tmp/candida_work:/tmp/candida_work "
       "quay.io/biocontainers/bcftools:1.24--h118bc1c_2 bcftools" % (os.getuid(), os.getgid()))
ref_nt = str(next(SeqIO.parse(REG, "fasta")).seq)
ref_prot = str(Seq(ref_nt).reverse_complement().translate()).rstrip("*")
def prot_from(vcf, hap):
    o = f"/tmp/candida_work/erg11_{os.path.basename(vcf)}.{hap}.fa"
    r = subprocess.run(f"{BCF} consensus -f {REG} -H {hap} {vcf} -o {o} 2>/dev/null",
                       shell=True, capture_output=True)
    if not os.path.exists(o): return None
    try:
        s = str(next(SeqIO.parse(o, "fasta")).seq)
    except Exception: return None
    finally:
        try: os.remove(o)
        except OSError: pass
    if len(s) != len(ref_nt): return None
    return str(Seq(s).reverse_complement().translate()).rstrip("*")
rows = []
vcfs = sorted(glob.glob(f"{C}/vcf/*.vcf.gz"))
for i, v in enumerate(vcfs, 1):
    acc = os.path.basename(v)[:-7]
    p1, p2 = prot_from(v, 1), prot_from(v, 2)
    if p1 is None and p2 is None:
        rows.append(dict(accession=acc, status="fail", substitutions="", n_subs="")); continue
    subs = {}
    for p, tag in ((p1, "a1"), (p2, "a2")):
        if p is None or len(p) != len(ref_prot): continue
        for j in range(len(ref_prot)):
            if p[j] != ref_prot[j] and p[j] not in "X*":
                subs.setdefault(f"{ref_prot[j]}{j+1}{p[j]}", set()).add(tag)
    zyg = {k: ("hom" if len(v_) == 2 else "het") for k, v_ in subs.items()}
    rows.append(dict(accession=acc, status="ok",
                     n_subs=len(subs),
                     substitutions=";".join(f"{k}({zyg[k]})" for k in sorted(subs, key=lambda x: int(x[1:-1]))),
                     Y132F=int("Y132F" in subs), S154F=int("S154F" in subs),
                     Y257H=int("Y257H" in subs), G464S=int("G464S" in subs)))
    if i % 50 == 0: print(f"  {i}/{len(vcfs)}", flush=True)
with open(f"{C}/erg11_genotypes.csv", "w", newline="") as fh:
    cols = ["accession","status","n_subs","substitutions","Y132F","S154F","Y257H","G464S"]
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
print("wrote", f"{C}/erg11_genotypes.csv", len(rows), "isolates")
