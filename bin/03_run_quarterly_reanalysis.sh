#!/usr/bin/env bash
set -euo pipefail

# Quarterly OMIM-driven genomic reanalysis runner.
# Edit config/reanalysis.config.sh first, then run:
#   bash bin/03_run_quarterly_reanalysis.sh config/my_run.config.sh

CONFIG=${1:-config/reanalysis.config.sh}
if [[ ! -f "$CONFIG" ]]; then
  echo "[ERROR] config file not found: $CONFIG" >&2
  exit 1
fi
source "$CONFIG"

mkdir -p "$OUTDIR" logs ref

if [[ "${DOWNLOAD_HPO_REFS:-0}" == "1" ]]; then
  echo "[STEP 0] Download/update HPO reference files"
  HPO_DOWNLOAD_ARGS=(--out-dir "${REF_DIR:-ref}" --hp-obo "${HP_OBO:-ref/hp.obo}" --genes-to-phenotype "${GENES_TO_PHENOTYPE:-ref/genes_to_phenotype.txt}")
  if [[ "${FORCE_DOWNLOAD_HPO:-0}" == "1" ]]; then HPO_DOWNLOAD_ARGS+=(--force); fi
  if [[ -n "${PHENOTYPE_TO_GENES:-}" ]]; then HPO_DOWNLOAD_ARGS+=(--phenotype-to-genes "$PHENOTYPE_TO_GENES"); fi
  python3 "$PIPELINE_DIR/bin/00_download_hpo_refs.py" "${HPO_DOWNLOAD_ARGS[@]}" 2> "$OUTDIR/00_download_hpo_refs.log"
fi

DATE_ARGS=()
if [[ -n "${START_MONTH:-}" ]]; then DATE_ARGS+=(--start "$START_MONTH"); fi
if [[ -n "${END_MONTH:-}" ]]; then DATE_ARGS+=(--end "$END_MONTH"); fi

UPDATES_TSV="$OUTDIR/omim_updates_raw.tsv"
GENE_TABLE="$OUTDIR/omim_new_gene_phenotype_table.tsv"
GENE_LIST="$OUTDIR/omim_new_gene_list.txt"

COLLECT_ARGS=(--out "$UPDATES_TSV" --sections "${OMIM_SECTIONS:-New Entries,New Clinical Synopses}")
if [[ "${USE_WEB:-1}" == "1" ]]; then COLLECT_ARGS+=(--web); fi
if [[ -n "${EMAIL_FILES:-}" ]]; then COLLECT_ARGS+=(--email-files ${EMAIL_FILES}); fi
if [[ "${USE_IMAP:-0}" == "1" ]]; then
  COLLECT_ARGS+=(--imap --imap-server "${IMAP_SERVER:-imap.163.com}" --imap-user "$IMAP_USER" --imap-password-env "${IMAP_PASSWORD_ENV:-OMIM_MAIL_PASSWORD}")
  if [[ -n "${IMAP_SINCE:-}" ]]; then COLLECT_ARGS+=(--imap-since "$IMAP_SINCE"); fi
fi

echo "[STEP 1] Collect OMIM updates"
python3 "$PIPELINE_DIR/bin/01_collect_omim_updates.py" \
  "${DATE_ARGS[@]}" \
  "${COLLECT_ARGS[@]}" 2> "$OUTDIR/01_collect.log"

if [[ ! -s "${GENES_TO_PHENOTYPE:-ref/genes_to_phenotype.txt}" ]]; then
  echo "[ERROR] genes_to_phenotype file not found: ${GENES_TO_PHENOTYPE:-ref/genes_to_phenotype.txt}" >&2
  echo "        Set DOWNLOAD_HPO_REFS=1 in the config or run bin/00_download_hpo_refs.py first." >&2
  exit 1
fi

echo "[STEP 2] Build OMIM gene-phenotype table"
python3 "$PIPELINE_DIR/bin/02_build_omim_gene_phenotype_table.py" \
  --updates "$UPDATES_TSV" \
  --genes-to-phenotype "${GENES_TO_PHENOTYPE:-ref/genes_to_phenotype.txt}" \
  --out "$GENE_TABLE" \
  --gene-list-out "$GENE_LIST" 2> "$OUTDIR/02_build_gene_table.log"

SNV_FOR_FILTER="${SNV_WITH_INHOUSE_AF:-}"

if [[ "${RUN_VARIANT_EXTRACTION:-0}" == "1" ]]; then
  echo "[STEP 3] Extract rare variants in updated OMIM genes"
  python3 "$EXTRACT_RARE_VARIANTS_SCRIPT" \
    --multianno "$MULTIANNO" \
    --sample-id "$REANALYSIS_SAMPLE_ID" \
    --gene-table "$GENE_TABLE" \
    --max-af "${GLOBAL_AF_CUTOFF:-0.01}" \
    --out-prefix "$OUTDIR/reanalysis_omim_new"

  RAW_EXTRACT="$OUTDIR/reanalysis_omim_new.rare_variants_in_OMIM2024_2025.tsv"
  SNV_FOR_FILTER="$RAW_EXTRACT"

  if [[ -n "${ADD_INHOUSE_AF_SCRIPT:-}" ]]; then
    echo "[STEP 4] Add in-house AF"
    python3 "$ADD_INHOUSE_AF_SCRIPT" \
      --input "$RAW_EXTRACT" \
      --sample-id "$REANALYSIS_SAMPLE_ID" \
      --output "$OUTDIR/reanalysis_omim_new.rare_variants.with_inhouseAF.tsv"
    SNV_FOR_FILTER="$OUTDIR/reanalysis_omim_new.rare_variants.with_inhouseAF.tsv"
  fi
fi

if [[ "${RUN_R_FILTER:-0}" == "1" ]]; then
  echo "[STEP 5] Filter candidate variants by inheritance and HPO matching"
  if [[ -z "$SNV_FOR_FILTER" ]]; then
    echo "[ERROR] RUN_R_FILTER=1 but SNV_WITH_INHOUSE_AF is empty and RUN_VARIANT_EXTRACTION=0" >&2
    exit 1
  fi
  Rscript "${R_FILTER_SCRIPT:-$PIPELINE_DIR/bin/04_filter_candidate_variants.R}" \
    --omim-gene-table "$GENE_TABLE" \
    --snv "$SNV_FOR_FILTER" \
    --genes-to-phenotype "${GENES_TO_PHENOTYPE:-ref/genes_to_phenotype.txt}" \
    --hp-obo "${HP_OBO:-hp.obo}" \
    --negative-id "$NEGATIVE_ID" \
    --sampleinfo "$SAMPLEINFO" \
    --pheno-xlsx "$PHENO_XLSX" \
    --out-prefix "$OUTDIR/filtered_variants" \
    --ar-freq-cutoff "${AR_FREQ_CUTOFF:-0.01}" \
    --ad-freq-cutoff "${AD_FREQ_CUTOFF:-0.0001}" \
    --global-af-cutoff "${GLOBAL_AF_CUTOFF:-0.01}" \
    --negative-denominator "${NEGATIVE_DENOMINATOR:-4128}" \
    --match-level "${HPO_MATCH_LEVEL:-exact}" \
    --min-hpo-intersect "${MIN_HPO_INTERSECT:-1}"
fi

echo "[DONE] OMIM reanalysis update finished. Output: $OUTDIR"
