# fungal-constraint-target-prioritization

Analysis code and derived data for:

> Vichathai S, Doungdun T, Ampawong S, Reamtong O, Huang Y-T, Sukphopetch P.
> *Evolutionary constraint reveals lineage-specific antifungal target prioritization and
> host-selectivity trade-offs across fungal pathogens.* Submitted to mSystems, 2026.

Archived on Zenodo: concept DOI [10.5281/zenodo.23386570](https://doi.org/10.5281/zenodo.23386570)
(always resolves to the latest release).

## Contents

| Folder | What it holds |
|---|---|
| `analysis_code/` | The 67 R, Python and shell scripts used for every result in the paper, at the exact version used. Names match the references in Text S1. |
| `environment/` | `renv.lock` (R 4.3.3 packages), `python_requirements.txt` (Python 3.12.3 packages), `micromamba_env_candida.yml` and `micromamba_env_gwas.yml` (read-mapping, variant-calling and pyseer environments), `R_sessionInfo.txt`, `tools_versions.md`. |
| `supplemental_data_files/` | Tables S1, S5, S9, S12–S17 and S19 as Excel workbooks, plus `data_dictionary.pdf`. These are the derived data cited in Data, Metadata and Code Availability. |

Inputs are public: NCBI assemblies (accession lists are produced by `download_ncbi_aspergillus.sh`
and `filter_assemblies_by_coverage.R`), SRA read sets for *C. tropicalis*
(`candida_process_isolate.sh`), the Af293 (GCF_000002655.1) and MYA-3404 (GCF_000006335.3)
references, UniProt reference proteomes, AlphaFold DB models, the PDB, KEGG and Gene Ontology
(retrieved 12 October 2026 and cached by the enrichment script), and the six published
*C. albicans* essentiality sources assembled by `build_calb_essentiality.py`.

## Script map (manuscript section → scripts)

| Step | Scripts |
|---|---|
| Genome retrieval and filtering | `download_ncbi_aspergillus.sh`, `filter_assemblies_by_coverage.R` |
| MIC phenotypes | `assemble_mic_datasheet.R`, `build_voriconazole_phenotype.R`, `build_mic_phenotype_numeric.R` |
| Association analysis (pyseer, 55 genomes) | `vrc_phase0_qc.R`, `vrc_phase2_snp_matrix.py`, `vrc_phase3_kinship.R`, `vrc_phase4_cyp51A_genotype.py`, `vrc_cyp51A_summary.R`, `vrc_phase5_report.R`, `vrc_phase6_candidate_panel.py`, `vrc_phase6_panel_test.R`, `vrc_phase6_confounding_check.R`, `vrc_panel_protein_pathway.R`, `vrc_panel_extended_genotype.py`, `vrc_weighting_build_genotypes.R`, `vrc_weighting_models.R`, `vrc_weighting_independence.R`, `vrc_team_analysis.R`, `vrc_case_series.R`, `vrc_cross_drug_pergene.R`, `vrc_promoter_repeat_scan.py` (TR34/TR46) |
| *A. fumigatus* variant extraction and πN/πS | `td_extract_variants.py`, `td_popstructure.py`, `td_constraint.py` |
| *C. tropicalis* processing and πN/πS | `candida_process_isolate.sh`, `candida_erg11_genotype.sh`, `candida_erg11_call.py`, `candida_asp_erg11_residue_map.py`, `candida_qc_and_cnv.R`, `candida_validate_cnv.R`, `candida_build_phenotype.R`, `candida_combine_results.R`, `candida_report_figures.R`, `ctrop_constraint.py` |
| Filtering funnel, structure and chemistry | `td_funnel.R`, `td_figures.R`, `rev_c1_canonical_filter_audit.R` |
| Orthology and *C. albicans* essentiality benchmark | `build_calb_essentiality.py`, `crossgenus_intersection.R`, `crossgenus_essentiality.R`, `rev_c0_calb_essentiality_per_locus.R`, `rev_b1_essentiality_benchmark.R` |
| Robustness and lineage subsampling | `rev_b2_asp_recompute.py`, `rev_b2_ctrop_recompute.py`, `rev_b2_analyse.R`, `rev_b2_crossgenus.R`, `rev_b2_figure.R` |
| Pathway enrichment | `rev_b3_pathway_enrichment.R` |
| Gene-family reconstruction and orthology grading | `rev_b4_01_assemble_family.sh`, `rev_b4_02_align_trim_tree.sh`, `rev_b4_02a_coverage.R`, `rev_b4_02b_constraints.R`, `rev_b4_02c_domain_tree.sh`, `rev_b4_03_tree_verdict_figure.R`, `rev_b4_04_synteny.R`, `rev_b4_05a_tier1_searches.sh`, `rev_b4_05_tier1_orthogroup_check.R`, `rev_b4_06_verdict_table.R` |
| Power simulation | `rev_b5_power_simulation.R` |
| Candidate-evidence table (Table S19) and main Table 1 | `rev_c2_00_candidate_set_union.R`, `rev_c2_candidate_evidence_master.R`, `rev_c3_main_table1.R` |
| Figures and supplemental data files | `ms_figures_tables.R`, `rev_fig5_benchmark_stability.R`, `build_supplemental_data_files.py` |

## Reproducing

1. Restore the R environment: `Rscript -e 'renv::restore(lockfile = "environment/renv.lock")'`.
   Install the Python packages with `pip install -r environment/python_requirements.txt`, and the
   two micromamba environments with `micromamba env create -f environment/micromamba_env_<name>.yml`.
   Command-line tool versions are listed in `environment/tools_versions.md`.
2. The scripts were run from a project root that is set near the top of each script
   (`B <- "/mnt/HuangLab_26T/HuangLab/Data/Tar_Nine/Ta"` in R, `BASE=` in shell, a `BASE`/`B`
   constant in Python). Edit that one line, or mount or symlink your checkout at that path. Scripts
   expect the layout `raw_data/`, `analysis/`, `analysis_code/`, `output/figure/`, `output/table/`.
3. Run the scripts in the order of the script map. Each writes CSV tables to `output/table/` and PDF
   figures to `output/figure/`. All random seeds are set inside the scripts.
4. Large intermediate outputs (alignments, VCFs, per-genome variant tables) are not in this
   repository; they are regenerated from the public inputs by the scripts in rows 1, 4 and 5.

## Licence

Code is released under the MIT License and the data files under CC BY 4.0 (see `LICENSE`).
