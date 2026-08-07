# OMIM update-driven genomic reanalysis pipeline

This repository provides a reproducible workflow for automated genomic reanalysis based on newly curated OMIM disease-gene associations. It is designed for quarterly or user-defined reanalysis of previously unresolved sequencing cases.

The pipeline generates **candidate variant review tables**. It does not replace manual clinical review, Sanger validation, segregation analysis, or ACMG/AMP-based interpretation.

## Workflow

1. Download HPO reference files: `hp.obo` and `genes_to_phenotype.txt`.
2. Collect OMIM updates from monthly OMIM update pages, saved notification emails, or an IMAP mailbox.
3. Extract disease MIM entries from `New Entries` and `New Clinical Synopses`.
4. Link disease MIM entries to gene symbols, inheritance modes, and HPO annotations using HPO `genes_to_phenotype.txt`.
5. Optionally extract rare variants in updated genes from initially unresolved cases.
6. Optionally annotate in-house carrier and allele frequencies.
7. Optionally filter candidate variants by inheritance model, allele frequency, consequence, gene constraint, damaging scores, and HPO phenotype matching.
8. Export inheritance-stratified candidate tables for manual review.

## Repository structure

```text
omim_reanalysis_pipeline_v3/
├── bin/
│   ├── 00_download_hpo_refs.py
│   ├── 01_collect_omim_updates.py
│   ├── 02_build_omim_gene_phenotype_table.py
│   ├── 03_run_quarterly_reanalysis.sh
│   └── 04_filter_candidate_variants.R
├── config/
│   └── reanalysis.config.example.sh
├── docs/
│   └── workflow.md
├── example/
├── ref/
├── results/
├── logs/
├── environment.yml
├── requirements.txt
├── .gitignore
└── README.md
```

## Installation

Python scripts use the Python standard library only. The R filtering script requires `dplyr`, `tidyr`, `stringr`, `readxl`, `writexl`, and `ontologyIndex`.

Using conda:

```bash
conda env create -f environment.yml
conda activate omim-reanalysis
```

Or install R packages manually:

```r
install.packages(c("dplyr", "tidyr", "stringr", "readxl", "writexl", "ontologyIndex"))
```

## Download HPO reference files

The pipeline does not bundle HPO reference files. Download the current release with:

```bash
python3 bin/00_download_hpo_refs.py --out-dir ref
```

This creates:

```text
ref/hp.obo
ref/genes_to_phenotype.txt
ref/download_manifest.tsv
```

To force a fresh download:

```bash
python3 bin/00_download_hpo_refs.py --out-dir ref --force
```

## Basic quarterly run

Copy and edit the example config:

```bash
cp config/reanalysis.config.example.sh config/my_run.config.sh
# edit PIPELINE_DIR, OUTDIR, reference paths, and downstream input paths
bash bin/03_run_quarterly_reanalysis.sh config/my_run.config.sh
```

Set a custom time range in the config:

```bash
START_MONTH=2026-01
END_MONTH=2026-03
```

Leave both empty to automatically use the previous completed quarter.

## Step 1: Collect OMIM updates from web

```bash
python3 bin/01_collect_omim_updates.py \
  --web \
  --start 2026-01 \
  --end 2026-03 \
  --sections "New Entries,New Clinical Synopses" \
  --out results/omim_updates_raw.tsv
```

## Step 1 alternative: saved email files

Save OMIM notification emails as `.txt` or `.eml`, then run:

```bash
python3 bin/01_collect_omim_updates.py \
  --email-files mail1.eml mail2.txt \
  --out results/omim_updates_from_email.tsv
```

## Step 1 alternative: IMAP mailbox

Use an app-specific mail authorization code through an environment variable. Do not write passwords into config files.

```bash
export OMIM_MAIL_PASSWORD='your-app-specific-password'

python3 bin/01_collect_omim_updates.py \
  --imap \
  --imap-server imap.example.com \
  --imap-user user@example.com \
  --imap-since 01-Jan-2026 \
  --out results/omim_updates_from_imap.tsv
```

## Step 2: Build updated OMIM gene-phenotype table

```bash
python3 bin/02_build_omim_gene_phenotype_table.py \
  --updates results/omim_updates_raw.tsv \
  --genes-to-phenotype ref/genes_to_phenotype.txt \
  --out results/omim_new_gene_phenotype_table.tsv \
  --gene-list-out results/omim_new_gene_list.txt
```

Main outputs:

```text
results/omim_new_gene_phenotype_table.tsv
results/omim_new_gene_list.txt
results/omim_new_gene_phenotype_table.tsv.unmatched_updates.tsv
```

## Step 3: Candidate variant filtering

After rare-variant extraction and in-house AF annotation, run:

```bash
Rscript bin/04_filter_candidate_variants.R \
  --omim-gene-table results/omim_new_gene_phenotype_table.tsv \
  --snv results/reanalysis_omim_new.rare_variants.with_inhouseAF.tsv \
  --genes-to-phenotype ref/genes_to_phenotype.txt \
  --hp-obo ref/hp.obo \
  --negative-id /path/to/unresolved_sample_ids.txt \
  --sampleinfo /path/to/sampleinfo_hpo.tsv \
  --pheno-xlsx /path/to/clinical_registry.xlsx \
  --out-prefix results/filtered_variants
```

Main outputs:

```text
results/filtered_variants.AD_XLD_candidates.tsv
results/filtered_variants.AR_XLR_candidates.tsv
results/filtered_variants.other_inheritance_candidates.tsv
results/filtered_variants.phenotype_matched_variants.tsv
results/filtered_variants.candidate_variants.xlsx
results/filtered_variants.qc_summary.tsv
```

## Configuration options

Important options in `config/reanalysis.config.example.sh`:

```bash
DOWNLOAD_HPO_REFS=1
FORCE_DOWNLOAD_HPO=0
GENES_TO_PHENOTYPE=ref/genes_to_phenotype.txt
HP_OBO=ref/hp.obo

RUN_VARIANT_EXTRACTION=0
RUN_R_FILTER=0

GLOBAL_AF_CUTOFF=0.01
AR_FREQ_CUTOFF=0.01
AD_FREQ_CUTOFF=0.0001
HPO_MATCH_LEVEL=exact
MIN_HPO_INTERSECT=1
```

## Data privacy

Patient-level clinical data, sample identifiers, raw VCF/BAM files, full cohort variant tables, mail passwords, and clinical registry files are not included in this repository. The `.gitignore` file is configured to reduce the chance of accidentally committing sensitive files, but users should still review all files before pushing to GitHub.

## Citation and reproducibility

For manuscript submission, create a GitHub release and cite the tagged version used in the study, for example:

```bash
git tag -a v1.0.0 -m "Version used for manuscript submission"
git push origin v1.0.0
```

Suggested statement:

> The OMIM update-driven genomic reanalysis pipeline is available at GitHub. The version used in this study is archived as release v1.0.0. Patient-level clinical data and raw variant files are not included because of privacy restrictions.
