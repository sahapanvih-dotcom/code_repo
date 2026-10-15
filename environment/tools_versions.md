# Command-line tools, as reported in Materials and Methods

| Tool | Version used |
|---|---|
| NCBI Datasets CLI | 18.36.0 |
| minimap2 | 2.26 |
| fastANI | 1.33 |
| Mash | 2.3 |
| BUSCO | 5 |
| MAFFT | 7.505 |
| trimAl | 1.5.rev1 |
| IQ-TREE | 2.0.7 |
| HMMER | 3.4 |
| BLAST+ | 2.12.0 |
| SeqKit | 2.13.0 |
| MUSCLE | 5.2 |
| pyseer | 1.4.0 (micromamba env `gwas`; the manuscript text states 1.3) |
| DIAMOND | 2.1.9 |
| fpocket | 4.0 |
| R | 4.3.3 (packages pinned in `renv.lock`) |
| Python | 3.12.3 (packages pinned in `python_requirements.txt`) |

Read mapping, variant calling and depth for *C. tropicalis* ran inside the micromamba environment
`candida` (`micromamba_env_candida.yml`); pyseer ran inside `gwas` (`micromamba_env_gwas.yml`).
Both YAML files are full exports of those environments.
