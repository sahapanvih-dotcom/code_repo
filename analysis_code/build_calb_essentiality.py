#!/usr/bin/env python3
"""Build a C. albicans SC5314 experimental gene-essentiality table from
pre-downloaded source files.  All inputs are passed as explicit paths.

Run with:  python3 -I build_calb_essentiality.py --out OUT.tsv --stats STATS.txt \
             --cgd ... --omeara ... --fu2021 ... --fu2025 ... --segal ...

No call is ever inferred from a gene name or from memory: every emitted row
carries the raw source string in the `evidence` column.
"""
import argparse
import collections
import csv
import os
import sys

# xlrd / openpyxl live in the user site-packages; -I drops that from sys.path.
# Append only that one trusted absolute path (never the data directory).
USER_SITE = "/home/huanglab/.local/lib/python3.12/site-packages"
if os.path.isdir(USER_SITE):
    sys.path.append(USER_SITE)

import openpyxl

OUT_COLS = ["gene_id", "gene_name", "call", "source", "evidence"]


def norm(v):
    return "" if v is None else str(v).strip()


def sheet_rows(path, sheet):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    rows = [[norm(c) for c in r] for r in ws.iter_rows(values_only=True)]
    wb.close()
    return rows


def get(row, i):
    return row[i] if i < len(row) else ""


# --------------------------------------------------------------------------
# 1. CGD phenotype_data.tab
# --------------------------------------------------------------------------
# Columns (1-based, per CGD download/phenotype/README):
#  1 Feature Name, 2 Feature Type, 3 Gene Name, 4 CGDID, 5 Reference,
#  6 Experiment Type, 7 Mutant Type, 8 Allele, 9 Strain Background,
# 10 Phenotype, ... 17 Species
CGD_KEEP_MUTANT = {"null", "repressible"}
CGD_MAP = {"inviable": "ESSENTIAL", "viable": "NONESSENTIAL"}


def parse_cgd(path, out, stats):
    kept = 0
    genes = set()
    skipped_mutant = collections.Counter()
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(rd)
        assert header[0] == "Feature Name" and header[9] == "Phenotype", header[:10]
        for row in rd:
            if len(row) < 10:
                continue
            pheno = row[9].strip()
            if pheno not in CGD_MAP:
                continue
            mutant = row[6].strip()
            if mutant not in CGD_KEEP_MUTANT:
                skipped_mutant[(mutant, pheno)] += 1
                continue
            gid = row[0].strip()
            if not gid:
                continue
            ev = ("Phenotype=%s | Mutant Type=%s | Experiment Type=%s | "
                  "Allele=%s | Strain Background=%s | Reference=%s | "
                  "Feature Type=%s" % (pheno, mutant, row[5].strip(),
                                       row[7].strip(), row[8].strip(),
                                       row[4].strip(), row[1].strip()))
            out.writerow([gid, row[2].strip(), CGD_MAP[pheno],
                          "CGD_phenotype", ev])
            kept += 1
            genes.add(gid)
    stats["CGD_phenotype"] = dict(rows=kept, genes=len(genes))
    stats["_cgd_skipped_mutant_types"] = dict(
        ("%s / %s" % k, v) for k, v in skipped_mutant.most_common())


# --------------------------------------------------------------------------
# 2. O'Meara 2015 (GRACE v1 tet-repression) ncomms7741-s2.xls
# --------------------------------------------------------------------------
OMEARA_MAP = {"E": "ESSENTIAL", "N": "NONESSENTIAL"}


def parse_omeara(path, out, stats):
    import xlrd
    bk = xlrd.open_workbook(path)
    sh = bk.sheet_by_name("Essentiality scores")
    hdr = [norm(c.value) for c in sh.row(4)]
    i_orf = hdr.index("orf19 name")
    i_com = hdr.index("Common")
    i_call = hdr.index("Essential/Not essential")
    i_dox = hdr.index("1 μg/mL DOX growth")
    i_tet = hdr.index("100 μg/mL tetracycline growth")
    i_foa = hdr.index("5-FOA excision")
    i_conc = hdr.index("Essentiality Concordance  (Y/N)")
    kept = 0
    genes = set()
    skipped = collections.Counter()
    for r in range(5, sh.nrows):
        row = [norm(c.value) for c in sh.row(r)]
        gid = get(row, i_orf)
        if not gid or gid == "orf19 name":
            continue
        raw = get(row, i_call)
        if raw not in OMEARA_MAP:
            skipped[raw] += 1
            continue
        ev = ("Essential/Not essential=%s | 1 ug/mL DOX growth=%s | "
              "100 ug/mL tetracycline growth=%s | 5-FOA excision=%s | "
              "Essentiality Concordance=%s" % (raw, get(row, i_dox),
                                               get(row, i_tet),
                                               get(row, i_foa),
                                               get(row, i_conc)))
        out.writerow([gid, get(row, i_com), OMEARA_MAP[raw],
                      "GRACEv1_OMeara2015", ev])
        kept += 1
        genes.add(gid)
    stats["GRACEv1_OMeara2015"] = dict(rows=kept, genes=len(genes),
                                       skipped_values=dict(skipped))


# --------------------------------------------------------------------------
# 3/4/5. GRACE verdict sheets (Fu 2021, Fu 2025) -- shared E/NE/GD vocabulary
# --------------------------------------------------------------------------
GRACE_MAP = {"E": "ESSENTIAL", "NE": "NONESSENTIAL", "GD": "UNKNOWN"}


def emit_grace(rows, hdr, cols, source, out, stats, verdict_idx=None,
               verdict_label="Essential_verdict", extra=""):
    """cols: dict with keys feature, orf, common, nodox, dox, verdict (header names)."""
    idx = {k: hdr.index(v) for k, v in cols.items()}
    if verdict_idx is not None:
        idx["verdict"] = verdict_idx
    kept = 0
    genes = set()
    skipped = collections.Counter()
    for row in rows:
        gid = get(row, idx["feature"])
        orf = get(row, idx["orf"])
        if not gid:
            gid = orf
        if not gid or gid == hdr[idx["feature"]]:
            continue
        raw = get(row, idx["verdict"])
        if raw not in GRACE_MAP:
            skipped[raw if raw else "<blank>"] += 1
            continue
        ev = ("%s=%s | No_DOX_score=%s | DOX_Score=%s | orf19=%s%s"
              % (verdict_label, raw, get(row, idx["nodox"]),
                 get(row, idx["dox"]), orf, extra))
        out.writerow([gid, get(row, idx["common"]), GRACE_MAP[raw],
                      source, ev])
        kept += 1
        genes.add(gid)
    stats[source] = dict(rows=kept, genes=len(genes),
                         skipped_values=dict(skipped))


def parse_fu2021(path, out, stats):
    # Sheet GRACE_Verdict  -> GRACE v1 rescored
    rows = sheet_rows(path, "GRACE_Verdict")
    hdr = rows[1]
    emit_grace(rows[2:], hdr,
               dict(feature="Feature Name", orf="orf19", common="Common",
                    nodox="No_DOX_score", dox="DOX_Score",
                    verdict="Essential_verdict"),
               "GRACEv1_Fu2021", out, stats,
               extra=" | screen=30 C YPD +/- DOX")
    # Sheet GRACEv2_Verdict
    rows = sheet_rows(path, "GRACEv2_Verdict")
    hdr = rows[1]
    emit_grace(rows[2:], hdr,
               dict(feature="Feature.Name", orf="ORF", common="Common.Name",
                    nodox="No_DOX_score", dox="DOX_Score",
                    verdict="Essential_verdict"),
               "GRACEv2_Fu2021", out, stats,
               extra=" | screen=30 C YPD +/- DOX")


def parse_fu2025(path, out, stats):
    """Primary screen scores: three stacked sections (GRACEv1/v2/v3).
    Only the GRACEv3 section is novel here (v1/v2 30 C columns are reprinted
    from Fu 2021, already captured from that primary source).
    30 C verdict is col 9; 39 C verdict is col 12."""
    rows = sheet_rows(path, "Primary screen scores")
    start = None
    for i, r in enumerate(rows):
        if get(r, 0) == "GRACEv3":
            start = i + 2          # skip section marker + column header row
            break
    if start is None:
        raise SystemExit("GRACEv3 section marker not found in Primary screen scores")
    hdr = rows[start - 1]
    assert get(hdr, 5) == "Feature.Name", hdr[:7]
    assert get(hdr, 9) == "Essential_verdict" and get(hdr, 12) == "Essential_verdict"
    emit_grace(rows[start:], hdr,
               dict(feature="Feature.Name", orf="ORF", common="Common.Name",
                    nodox="No_DOX_score", dox="DOX_Score"),
               "GRACEv3_Fu2025", out, stats, verdict_idx=9,
               extra=" | screen=30 C YPD +/- DOX (GRACEv3 section)")


# --------------------------------------------------------------------------
# 6. Segal 2018 transposon essentiality
# --------------------------------------------------------------------------
SEGAL_MAP = {"Ess": "ESSENTIAL", "NE": "NONESSENTIAL"}


def parse_segal(path, out, stats):
    rows = sheet_rows(path, "2A-Ca classification")
    hdr = rows[5]
    i_id = hdr.index("Standard name")
    i_com = hdr.index("Common name")
    i_call = hdr.index("CaTn")
    i_rf = hdr.index("Ca RF verdict CaTn")
    i_hits = hdr.index("Hits")
    i_reads = hdr.index("Reads")
    i_fi = hdr.index("Freedom index")
    i_ni = hdr.index("Neighborhood index")
    kept = 0
    genes = set()
    skipped = collections.Counter()
    for row in rows[6:]:
        gid = get(row, i_id)
        if not gid or gid == "Standard name":
            continue
        raw = get(row, i_call)
        if raw not in SEGAL_MAP:
            skipped[raw if raw else "<blank>"] += 1
            continue
        ev = ("CaTn=%s | Ca RF verdict CaTn (score)=%s | Hits=%s | Reads=%s | "
              "Freedom index=%s | Neighborhood index=%s"
              % (raw, get(row, i_rf), get(row, i_hits), get(row, i_reads),
                 get(row, i_fi), get(row, i_ni)))
        out.writerow([gid, get(row, i_com), SEGAL_MAP[raw], "Segal2018", ev])
        kept += 1
        genes.add(gid)
    stats["Segal2018"] = dict(rows=kept, genes=len(genes),
                              skipped_values=dict(skipped))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--cgd")
    ap.add_argument("--omeara")
    ap.add_argument("--fu2021")
    ap.add_argument("--fu2025")
    ap.add_argument("--segal")
    a = ap.parse_args()

    stats = {}
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        out = csv.writer(fh, delimiter="\t", lineterminator="\n",
                         quoting=csv.QUOTE_MINIMAL)
        out.writerow(OUT_COLS)
        if a.cgd:
            parse_cgd(a.cgd, out, stats)
        if a.omeara:
            parse_omeara(a.omeara, out, stats)
        if a.fu2021:
            parse_fu2021(a.fu2021, out, stats)
        if a.fu2025:
            parse_fu2025(a.fu2025, out, stats)
        if a.segal:
            parse_segal(a.segal, out, stats)

    with open(a.stats, "w", encoding="utf-8") as fh:
        for k, v in stats.items():
            fh.write("%s\t%s\n" % (k, v))
    for k, v in stats.items():
        print("%s\t%s" % (k, v))


if __name__ == "__main__":
    main()
