# Example configuration for OMIM update-driven reanalysis.
# Copy this file to config/my_run.config.sh and edit paths before running.

# Core paths
PIPELINE_DIR=/path/to/omim_reanalysis_pipeline_v3
OUTDIR=results/omim_reanalysis_$(date +%Y%m%d)

# Time range. Leave empty to use the previous completed quarter.
# Format: YYYY-MM, for example START_MONTH=2026-01 and END_MONTH=2026-03
START_MONTH=
END_MONTH=

# OMIM update source
USE_WEB=1
OMIM_SECTIONS="New Entries,New Clinical Synopses"

# Optional saved email text/eml files, space-separated.
EMAIL_FILES=

# Optional IMAP mode. Use an app-specific authorization code through an environment variable.
USE_IMAP=0
IMAP_SERVER=imap.example.com
IMAP_USER=user@example.com
IMAP_PASSWORD_ENV=OMIM_MAIL_PASSWORD
IMAP_SINCE=

# HPO reference files
REF_DIR=ref
DOWNLOAD_HPO_REFS=1
FORCE_DOWNLOAD_HPO=0
GENES_TO_PHENOTYPE=ref/genes_to_phenotype.txt
HP_OBO=ref/hp.obo
# Optional; leave empty unless needed.
PHENOTYPE_TO_GENES=

# Optional downstream variant extraction
RUN_VARIANT_EXTRACTION=0
EXTRACT_RARE_VARIANTS_SCRIPT=/path/to/01_extract_omim2024_2025_rare_variants.safe.py
ADD_INHOUSE_AF_SCRIPT=/path/to/02_add_inhouse_af.py
MULTIANNO=/path/to/cohort.hg19_multianno.txt
REANALYSIS_SAMPLE_ID=/path/to/unresolved_sample_ids.txt
GLOBAL_AF_CUTOFF=0.01

# Candidate filtering
RUN_R_FILTER=0
R_FILTER_SCRIPT=/path/to/omim_reanalysis_pipeline_v3/bin/04_filter_candidate_variants.R
# If RUN_VARIANT_EXTRACTION=0, provide a precomputed SNV table here.
SNV_WITH_INHOUSE_AF=
NEGATIVE_ID=/path/to/unresolved_sample_ids.txt
SAMPLEINFO=/path/to/sampleinfo_hpo.tsv
PHENO_XLSX=/path/to/clinical_registry.xlsx

AR_FREQ_CUTOFF=0.01
AD_FREQ_CUTOFF=0.0001
NEGATIVE_DENOMINATOR=4128
HPO_MATCH_LEVEL=exact
MIN_HPO_INTERSECT=1
