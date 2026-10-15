from Bio import SeqIO, Align
import re, subprocess
# A. fumigatus cyp51A
asp = {r.id: r for r in SeqIO.parse("doc/Aspergillus/voriconazole_analysis_2026-08/panel_annotation/panel_proteins.faa","fasta")}
cyp = [r for k,r in asp.items() if "752137" in k or "cyp51A" in r.description]
print("cyp51A candidates:", [ (r.id, len(r.seq)) for r in cyp ][:5])
# C. tropicalis ERG11 protein id from gff
gff = open("analysis/candida/ref/genomic.gff").read()
pid = None
for line in gff.split("\n"):
    if "CTRG_05283" in line and "\tCDS\t" in line:
        m = re.search(r"protein_id=([^;]+)", line)
        if m: pid = m.group(1); break
print("C. tropicalis ERG11 protein_id:", pid)
ctr = {r.id: r for r in SeqIO.parse("analysis/candida/ref/protein.faa","fasta")}
erg = ctr.get(pid)
print("ERG11 length:", len(erg.seq) if erg else None)
a = str(cyp[0].seq).rstrip("*"); b = str(erg.seq).rstrip("*")
print("cyp51A length:", len(a))
al = Align.PairwiseAligner(open_gap_score=-11, extend_gap_score=-1, substitution_matrix=Align.substitution_matrices.load("BLOSUM62"), mode="global")
aln = al.align(a,b)[0]
print("score:", aln.score)
# build position maps from aligned blocks
A,B = aln.aligned
asp2can, can2asp = {}, {}
for (as_,ae),(bs,be) in zip(A,B):
    for i in range(ae-as_):
        asp2can[as_+i+1] = bs+i+1
        can2asp[bs+i+1] = as_+i+1
ident = sum(1 for (as_,ae),(bs,be) in zip(A,B) for i in range(ae-as_) if a[as_+i]==b[bs+i])
print("aligned identity: %.1f%% over %d aligned residues" % (100*ident/sum(ae-as_ for (as_,ae),_ in zip(A,B)), sum(ae-as_ for (as_,ae),_ in zip(A,B))))
print("\n=== Candida ERG11 hotspot -> A. fumigatus cyp51A position")
for p,lab in [(132,"Y132F"),(143,"K143R"),(145,"F145L"),(154,"S154F"),(405,"S405F"),(446,"D446E"),(448,"G448E"),(449,"F449V"),(450,"G450E"),(464,"G464S")]:
    q = can2asp.get(p)
    print(f"  Ct {lab:7s} (res {b[p-1]}{p})  ->  Af cyp51A  {a[q-1]+str(q) if q else 'gap/unaligned'}")
print("\n=== A. fumigatus cyp51A hotspot -> Candida ERG11 position")
for p,lab in [(98,"L98H"),(121,"Y121F"),(289,"T289A"),(448,"G448S"),(54,"G54"),(138,"G138"),(220,"M220")]:
    q = asp2can.get(p)
    print(f"  Af {lab:6s} (res {a[p-1]}{p})  ->  Ct ERG11  {b[q-1]+str(q) if q else 'gap/unaligned'}")

# emit the verified mapping as a table
import csv as _csv
rows=[]
for pos,lab in [(132,"Y132F"),(143,"K143R"),(145,"F145L"),(154,"S154F"),(405,"S405F"),
                (446,"D446E"),(448,"G448E"),(449,"F449V"),(450,"G450E"),(464,"G464S")]:
    q=can2asp.get(pos)
    rows.append(dict(direction="Ct->Af", candida_allele=lab, candida_res=b[pos-1]+str(pos),
                     asp_res=(a[q-1]+str(q)) if q else "unaligned",
                     same_residue=bool(q and a[q-1]==b[pos-1])))
for pos,lab in [(98,"L98H"),(121,"Y121F"),(289,"T289A"),(448,"G448S"),(54,"G54"),(138,"G138"),(220,"M220")]:
    q=asp2can.get(pos)
    rows.append(dict(direction="Af->Ct", candida_allele=lab, candida_res=(b[q-1]+str(q)) if q else "unaligned",
                     asp_res=a[pos-1]+str(pos),
                     same_residue=bool(q and a[pos-1]==b[q-1])))
with open("doc/Candida/tables/erg11_cyp51A_residue_map.csv","w",newline="") as fh:
    w=_csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("\nwrote doc/Candida/tables/erg11_cyp51A_residue_map.csv", len(rows), "rows")
print("identity_pct", round(100*ident/sum(ae-as_ for (as_,ae),_ in zip(A,B)),1))
