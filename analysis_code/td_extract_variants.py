#!/usr/bin/env python3
"""Align one Aspergillus assembly to Af293 and extract substitutions + covered intervals.

Parses minimap2's cs tag directly (no paftools/k8 on this host). The cs string is
expressed in REFERENCE orientation, so walking it from the alignment's ref start
yields reference coordinates without strand bookkeeping.

Outputs, per genome, into analysis/target_discovery/var/:
  <acc>.snp.tsv.gz   contig, 0-based ref pos, ref base, alt base
  <acc>.cov.bed      reference intervals covered by a primary alignment
Coverage is tracked because absence of variants in an unaligned region is NOT
evidence of conservation - it is absence of data.
"""
import subprocess, sys, os, re, gzip

REF = "analysis/target_discovery/ref/Af293.fna"
OUT = "analysis/target_discovery/var"
MIN_MAPQ = 30
MIN_BLOCK = 2000          # ignore short spurious blocks

cs_re = re.compile(rb"([:*+-])([A-Za-z0-9]+)")

def run(genome_fa, acc):
    snp_path = os.path.join(OUT, acc + ".snp.tsv.gz")
    cov_path = os.path.join(OUT, acc + ".cov.bed")
    if os.path.exists(snp_path) and os.path.exists(cov_path):
        return "skip"
    cmd = ["minimap2", "-cx", "asm5", "--cs", "-t", "2", "--secondary=no", REF, genome_fa]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    snps, cov = [], []
    for line in proc.stdout:
        f = line.rstrip(b"\n").split(b"\t")
        if len(f) < 12:
            continue
        try:
            mapq = int(f[11])
        except ValueError:
            continue
        ctg, rs, re_ = f[5].decode(), int(f[7]), int(f[8])
        if mapq < MIN_MAPQ or (re_ - rs) < MIN_BLOCK:
            continue
        tags = {t[:5]: t for t in f[12:]}
        cs = next((t for t in f[12:] if t.startswith(b"cs:Z:")), None)
        if cs is None:
            continue
        cov.append((ctg, rs, re_))
        pos = rs
        for op, val in cs_re.findall(cs[5:]):
            if op == b":":                      # identical block
                pos += int(val)
            elif op == b"*":                    # substitution: ref then alt
                snps.append((ctg, pos, val[:1].decode().upper(), val[1:2].decode().upper()))
                pos += 1
            elif op == b"-":                    # deletion from reference
                pos += len(val)
            elif op == b"+":                    # insertion into reference
                pass
    proc.wait()
    with gzip.open(snp_path, "wt") as fh:
        for c, p, r, a in snps:
            fh.write(f"{c}\t{p}\t{r}\t{a}\n")
    with open(cov_path, "w") as fh:
        for c, s, e in sorted(cov):
            fh.write(f"{c}\t{s}\t{e}\n")
    return f"{len(snps)} snps, {len(cov)} blocks"

if __name__ == "__main__":
    genome_fa, acc = sys.argv[1], sys.argv[2]
    try:
        print(f"{acc}\t{run(genome_fa, acc)}", flush=True)
    except Exception as ex:
        print(f"{acc}\tFAIL {ex}", flush=True)
