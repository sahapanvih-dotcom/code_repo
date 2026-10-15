#!/usr/bin/env python3
"""Assemble the nine supplemental data files promised in Text S1 §5 of the mSystems submission.

Each Table S* bundles several generated CSV/TSV files, so each is written as one multi-sheet
.xlsx workbook (an ASM-accepted supplemental format). Python rather than R because the job is
file packaging, not analysis: pandas/openpyxl write multi-sheet workbooks in one pass.
The underlying CSVs in output/table/ remain the canonical lab-convention outputs.
"""
import os, re, gzip, sys
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(BASE, "doc/manuscript_2026-10/submission_mSystems/supplemental_data_files")
os.makedirs(OUT, exist_ok=True)
T = lambda *p: os.path.join(BASE, *p)

def rev(prefix):
    d = T("output/table")
    return [os.path.join(d, f) for f in sorted(os.listdir(d))
            if f.startswith(prefix + "_") and f.endswith(".csv")]

TABLES = {
 "Table_S1":  ("Assembly inventory and quality control",
               [T("output/table/assembly_coverage_filter.csv")]),
 "Table_S5":  ("Per-gene constraint, the canonical filter audit and the candidate lists",
               [T("analysis/target_discovery/constraint_per_gene_dedup.tsv"),
                T("analysis/candida/constraint/ctrop_constraint_per_gene.tsv"),
                T("output/table/rev_c1_canonical_filter_audit.csv"),
                T("output/table/rev_c1_funnel_counts_canonical.csv"),
                T("output/table/rev_c1_candidate_815_audit.csv"),
                T("output/table/ms_table7_crossgenus_candidates.csv")]),
 "Table_S9":  ("C. albicans essentiality calls and the cross-genus orthology map",
               [T("analysis/crossgenus/calb_essentiality.tsv.gz"),
                T("analysis/crossgenus/ortho_map.tsv"),
                T("analysis/crossgenus/calb_idmap.tsv")]),
 "Table_S12": ("Essentiality-benchmark definitions and per-model performance",
               [T("output/table/rev_b1_outcome_definitions.csv"),
                T("output/table/rev_b1_performance.csv")]),
 "Table_S13": ("Paired model differences and gene-family blocks",
               [T("output/table/rev_b1_paired_deltas.csv"),
                T("output/table/rev_b1_gene_family_blocks.csv")]),
 "Table_S14": ("Ranking robustness and the matched cross-genus comparison", rev("rev_b2")),
 "Table_S15": ("Pathway and host-homology enrichment", rev("rev_b3")),
 "Table_S16": ("Gene-family reconstruction and Tier-1 orthology confidence", rev("rev_b4")),
 "Table_S17": ("Association power simulation", rev("rev_b5")),
 "Table_S19": ("Candidate-level evidence and triage for all 26 lead candidates",
               [T("output/table/rev_c2_candidate_evidence_master.csv"),
                T("output/table/rev_c3_main_table1.csv"),
                T("output/table/rev_c1_tier1_reevaluation.csv"),
                T("output/table/rev_b4_tier1_orthology_confidence_union26.csv"),
                T("output/table/rev_b4_tier1_paralog_detail_union26.csv"),
                T("output/table/rev_c0_calb_essentiality_per_af_locus.csv"),
                T("output/table/rev_c0_essentiality_crosscheck.csv")]),
}

def read_any(p):
    op = gzip.open if p.endswith(".gz") else open
    sep = "\t" if re.search(r"\.tsv(\.gz)?$", p) else ","
    with op(p, "rt") as fh:
        return pd.read_csv(fh, sep=sep, dtype=str, keep_default_na=False)

def sheet_name(path, used):
    """Excel caps sheet names at 31 characters. Truncating from the right alone made the
    C. tropicalis sensitivity arm indistinguishable from the primary run, so the sensitivity
    suffix is moved to a short prefix and over-long names are elided in the middle, which
    keeps both ends readable. Every workbook also carries an _INDEX sheet mapping sheet to
    source file, so nothing depends on reading the abbreviation."""
    n = re.sub(r"\.(csv|tsv)(\.gz)?$", "", os.path.basename(path))
    n = re.sub(r"^rev_[bc]\d_", "", n)
    pre = ""
    if n.endswith("_ctshrunk_sensitivity"):
        n = n[: -len("_ctshrunk_sensitivity")]; pre = "SENS_"
    room = 31 - len(pre)
    if len(n) > room:
        head = (room - 1) * 2 // 3
        n = n[:head] + "~" + n[len(n) - (room - 1 - head):]
    n = pre + n
    i = 1
    while n in used:
        i += 1; n = (n[: 31 - len(str(i)) - 1] + "_%d" % i)
    used.add(n); return n

for name, (title, paths) in TABLES.items():
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        sys.exit("missing source for %s: %s" % (name, missing))
    dest = os.path.join(OUT, name + ".xlsx")
    used, rows = {"_INDEX"}, []
    frames = []
    for p in paths:
        df = read_any(p)
        s = sheet_name(p, used)
        frames.append((s, df))
        rows.append((s, os.path.relpath(p, BASE), len(df), len(df.columns)))
    index = pd.DataFrame([{"sheet": s, "source_file": src, "rows": nr, "columns": nc}
                          for s, src, nr, nc in rows])
    index.insert(0, "table", name)
    index.insert(1, "contents", title)
    with pd.ExcelWriter(dest, engine="openpyxl") as xw:
        index.to_excel(xw, sheet_name="_INDEX", index=False)
        for s, df in frames:
            df.to_excel(xw, sheet_name=s, index=False)
    print("%-10s %-55s %2d sheets %8d rows" % (name, title, len(rows), sum(r[2] for r in rows)))
    for s, src, nr, nc in rows:
        print("             %-32s <- %-60s %7d x %3d" % (s, src, nr, nc))
