#!/usr/bin/env Rscript

# OMIM-driven candidate SNV filtering
#
# This script is a parameterized, GitHub-ready version of the original
# R_find_candidatesnv2.R workflow. It takes an OMIM updated gene table,
# an extracted SNV table, HPO gene-to-phenotype annotations, patient HPO
# annotations, and optionally a clinical registry Excel file, then exports
# AD/XLD, AR/XLR, and other-inheritance candidate variant tables.

suppressPackageStartupMessages({
  required_pkgs <- c("dplyr", "tidyr", "stringr", "ontologyIndex", "readxl")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))",
      call. = FALSE
    )
  }
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ontologyIndex)
  library(readxl)
})

# -----------------------------
# Minimal command-line parser
# -----------------------------
parse_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  opts <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    key <- sub("^--", "", key)
    key <- gsub("-", "_", key)
    if (i == length(argv) || startsWith(argv[[i + 1]], "--")) {
      opts[[key]] <- TRUE
      i <- i + 1
    } else {
      opts[[key]] <- argv[[i + 1]]
      i <- i + 2
    }
  }
  opts
}

usage <- function() {
  cat("\nUsage:\n")
  cat("  Rscript bin/04_filter_candidate_variants.R \\\n")
  cat("    --omim-gene-table results/omim_new_gene_phenotype_table.tsv \\\n")
  cat("    --snv results/reanalysis_omim_new.rare_variants.with_inhouseAF.tsv \\\n")
  cat("    --genes-to-phenotype ref/genes_to_phenotype.txt \\\n")
  cat("    --hp-obo ref/hp.obo \\\n")
  cat("    --negative-id reanalysis.id \\\n")
  cat("    --sampleinfo sampleinfo_diagnosis_level2_0411.tsv \\\n")
  cat("    --out-prefix results/filtered_variants\n\n")
}

opts <- parse_args()
if (isTRUE(opts$help) || isTRUE(opts$h)) {
  usage()
  quit(save = "no", status = 0)
}

get_opt <- function(name, default = NULL, required = FALSE) {
  key <- gsub("-", "_", name)
  val <- opts[[key]]
  if (is.null(val) || identical(val, "")) {
    if (required) stop("Missing required option --", gsub("_", "-", key), call. = FALSE)
    return(default)
  }
  val
}

# -----------------------------
# Options
# -----------------------------
omim_gene_table <- get_opt("omim-gene-table", required = TRUE)
snv_file <- get_opt("snv", required = TRUE)
genes_to_phenotype_file <- get_opt("genes-to-phenotype", required = TRUE)
hp_obo_file <- get_opt("hp-obo", required = TRUE)
negative_id_file <- get_opt("negative-id", default = NULL)
sampleinfo_file <- get_opt("sampleinfo", required = TRUE)
pheno_xlsx <- get_opt("pheno-xlsx", default = NULL)
out_prefix <- get_opt("out-prefix", required = TRUE)

clinical_sample_col <- get_opt("clinical-sample-col", default = "样本号")
sample_id_col <- get_opt("sample-id-col", default = "SampleID")

negative_freq_col_opt <- get_opt("negative-freq-col", default = "auto")
positive_freq_col_opt <- get_opt("positive-freq-col", default = "auto")
global_af_col <- get_opt("global-af-col", default = "AF")

ar_freq_cutoff <- as.numeric(get_opt("ar-freq-cutoff", default = "0.01"))
ad_freq_cutoff <- as.numeric(get_opt("ad-freq-cutoff", default = "0.0001"))
global_af_cutoff <- as.numeric(get_opt("global-af-cutoff", default = as.character(ar_freq_cutoff)))
revel_cutoff <- as.numeric(get_opt("revel-cutoff", default = "0.644"))
cadd_cutoff <- as.numeric(get_opt("cadd-cutoff", default = "25"))
pli_cutoff <- as.numeric(get_opt("pli-cutoff", default = "0.9"))
negative_denominator <- as.numeric(get_opt("negative-denominator", default = "4128"))
min_hpo_intersect <- as.integer(get_opt("min-hpo-intersect", default = "1"))
match_level <- get_opt("match-level", default = "exact")

if (!match_level %in% c("exact", "expanded")) {
  stop("--match-level must be either 'exact' or 'expanded'", call. = FALSE)
}

out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------
# Helper functions
# -----------------------------
message2 <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")

safe_read_tsv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_squish(x)
  x <- str_replace(x, "\\.0$", "")
  x
}

to_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  x <- as.character(x)
  x[x %in% c("", ".", "NA", "NaN", "nan", "null", "NULL")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

first_existing <- function(df, candidates, required = FALSE, what = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[[1]])
  if (required) stop("Cannot find ", what, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  NULL
}

rename_if_present <- function(df, old, new) {
  if (old %in% colnames(df) && !(new %in% colnames(df))) {
    colnames(df)[colnames(df) == old] <- new
  }
  df
}

split_ids <- function(x) {
  if (is.na(x) || x == "" || x == ".") return(character())
  ids <- unlist(str_split(as.character(x), "[,;]"))
  ids <- str_trim(ids)
  ids <- ids[ids != "" & ids != "."]
  ids[ids == "HP:0025356"] <- "HP:0001263"  # obsolete psychomotor retardation -> GDD fallback
  unique(ids)
}

collapse_unique <- function(x, sep = ";") {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "."]
  if (length(x) == 0) return(".")
  paste(sort(unique(x)), collapse = sep)
}

intersect_count <- function(a, b) {
  aa <- split_ids(a)
  bb <- split_ids(b)
  length(intersect(aa, bb))
}

intersect_list <- function(a, b) {
  aa <- split_ids(a)
  bb <- split_ids(b)
  hit <- intersect(aa, bb)
  if (length(hit) == 0) return("")
  paste(hit, collapse = ",")
}

hpo_ids_to_terms <- function(x, hpo_name_map) {
  ids <- split_ids(x)
  if (length(ids) == 0) return(NA_character_)
  terms <- vapply(ids, function(id) {
    if (id %in% names(hpo_name_map)) hpo_name_map[[id]] else paste0(id, " (unknown)")
  }, character(1))
  paste(terms, collapse = ", ")
}

parse_gt <- function(gt_or_info) {
  if (is.na(gt_or_info) || gt_or_info == "" || gt_or_info == ".") return(".")
  strsplit(as.character(gt_or_info), ":", fixed = TRUE)[[1]][[1]]
}

is_hom_or_hemi_alt <- function(gt) {
  gt <- parse_gt(gt)
  if (gt == ".") return(FALSE)
  gt2 <- gsub("\\|", "/", gt)
  alleles <- unlist(strsplit(gt2, "/", fixed = TRUE))
  alleles <- alleles[alleles != "." & alleles != ""]
  if (length(alleles) == 0) return(FALSE)
  all(alleles != "0")
}

has_inheritance <- function(x, codes) {
  vals <- unlist(strsplit(as.character(x), "[,;|/ ]+"))
  vals <- vals[vals != ""]
  any(vals %in% codes)
}

is_lof_variant <- function(func, exonic_func, lof_col = NULL) {
  func <- tolower(as.character(func))
  exonic_func <- tolower(as.character(exonic_func))
  lof_flag <- rep(FALSE, length(func))
  if (!is.null(lof_col)) {
    lof_flag <- !is.na(lof_col) & lof_col != "" & lof_col != "."
  }
  lof_func <- exonic_func %in% c(
    "frameshift deletion", "frameshift insertion", "frameshift substitution",
    "stopgain", "stoploss", "startloss"
  )
  splicing <- func == "splicing" | grepl("splicing", func)
  lof_flag | lof_func | splicing
}

is_protein_altering_variant <- function(func, exonic_func) {
  func <- tolower(as.character(func))
  exonic_func <- tolower(as.character(exonic_func))
  splicing <- func == "splicing" | grepl("splicing", func)
  ncrna <- func == "ncrna_exonic"
  exonic_altering <- grepl("exonic", func) & !(exonic_func %in% c("synonymous snv", "unknown", ".", ""))
  splicing | ncrna | exonic_altering
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_excel_if_possible <- function(sheet_list, path) {
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheet_list, path = path)
    return(TRUE)
  }
  warning("Package 'writexl' not installed. Excel output skipped: ", path)
  FALSE
}

# -----------------------------
# Load HPO ontology
# -----------------------------
message2("Reading HPO ontology: ", hp_obo_file)
hpo <- ontologyIndex::get_ontology(hp_obo_file, extract_tags = "everything")
hpo_name_map <- hpo$name
root_term <- "HP:0000118"
level2_terms <- names(hpo$parents)[vapply(hpo$parents, function(parents) root_term %in% parents, logical(1))]

get_level2_ancestors <- function(term, ontology, level2_set, visited = character()) {
  if (is.na(term) || term == "" || !(term %in% names(ontology$parents))) return(character())
  if (term %in% visited) return(character())
  visited <- c(visited, term)
  parents <- ontology$parents[[term]]
  parents <- parents[!is.na(parents)]
  result <- character()
  for (p in parents) {
    if (p %in% level2_set) {
      result <- c(result, p)
    } else {
      result <- c(result, get_level2_ancestors(p, ontology, level2_set, visited))
    }
  }
  unique(result)
}

# -----------------------------
# OMIM updated gene table
# -----------------------------
message2("Reading OMIM updated gene table: ", omim_gene_table)
omim <- safe_read_tsv(omim_gene_table)
omim <- rename_if_present(omim, "Gene.Locus", "Gene/Locus")
omim <- rename_if_present(omim, "Gene", "Gene/Locus")
omim <- rename_if_present(omim, "gene_symbol", "Gene/Locus")

gene_col <- first_existing(omim, c("Gene/Locus", "gene_id", "Matched_OMIM_gene"), required = TRUE, what = "OMIM gene column")
phenotype_col <- first_existing(omim, c("Phenotype", "OMIM_Phenotype", "Title"), required = FALSE)
inheritance_col <- first_existing(omim, c("Inheritance", "OMIM_Inheritance"), required = FALSE)
mim_col <- first_existing(omim, c("MIM_number", "MIM", "MIM_number.x", "MIM Number"), required = FALSE)
update_type_col <- first_existing(omim, c("Update_type", "Section"), required = FALSE)
update_time_col <- first_existing(omim, c("update time", "Date", "OMIM_update_time"), required = FALSE)

omim_core <- omim %>%
  mutate(
    gene_id = clean_id(.data[[gene_col]]),
    OMIM_Phenotype = if (!is.null(phenotype_col)) as.character(.data[[phenotype_col]]) else ".",
    OMIM_Inheritance = if (!is.null(inheritance_col)) as.character(.data[[inheritance_col]]) else ".",
    OMIM_MIM_number = if (!is.null(mim_col)) as.character(.data[[mim_col]]) else ".",
    OMIM_Update_type = if (!is.null(update_type_col)) as.character(.data[[update_type_col]]) else ".",
    OMIM_update_time = if (!is.null(update_time_col)) as.character(.data[[update_time_col]]) else "."
  ) %>%
  filter(!is.na(gene_id), gene_id != "", gene_id != ".") %>%
  group_by(gene_id) %>%
  summarise(
    OMIM_Phenotype = collapse_unique(OMIM_Phenotype),
    OMIM_Inheritance = collapse_unique(OMIM_Inheritance),
    OMIM_MIM_number = collapse_unique(OMIM_MIM_number),
    OMIM_Update_type = collapse_unique(OMIM_Update_type),
    OMIM_update_time = collapse_unique(OMIM_update_time),
    .groups = "drop"
  )

target_genes <- sort(unique(omim_core$gene_id))
target_mims <- unique(gsub("[^0-9]", "", unlist(strsplit(paste(omim_core$OMIM_MIM_number, collapse = ";"), "[;,]"))))
target_mims <- target_mims[target_mims != ""]
message2("Target OMIM genes: ", length(target_genes))

# -----------------------------
# HPO genes_to_phenotype
# -----------------------------
message2("Reading genes_to_phenotype: ", genes_to_phenotype_file)
g2p_raw <- safe_read_tsv(genes_to_phenotype_file)
col_lower <- tolower(colnames(g2p_raw))
colnames(g2p_raw) <- make.names(colnames(g2p_raw), unique = TRUE)
# Use original/lower lookup after make.names as well
lookup <- setNames(colnames(g2p_raw), tolower(gsub("\\.", "_", colnames(g2p_raw))))
find_g2p_col <- function(candidates, fallback_index = NULL) {
  norm_candidates <- tolower(gsub("[- .]", "_", candidates))
  for (cand in norm_candidates) {
    if (cand %in% names(lookup)) return(lookup[[cand]])
  }
  if (!is.null(fallback_index) && ncol(g2p_raw) >= fallback_index) return(colnames(g2p_raw)[fallback_index])
  stop("Cannot detect genes_to_phenotype column. Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
}

g2p_gene_col <- find_g2p_col(c("gene_symbol", "gene", "entrez_gene_symbol", "gene_id"), 2)
g2p_hpo_id_col <- find_g2p_col(c("hpo_id", "hpo"), 3)
g2p_hpo_name_col <- find_g2p_col(c("hpo_name", "term_name"), 4)
g2p_freq_col <- find_g2p_col(c("frequency", "freq"), 5)
g2p_disease_col <- find_g2p_col(c("disease_id", "disease"), 6)

gene2hpo <- g2p_raw %>%
  transmute(
    gene_id = clean_id(.data[[g2p_gene_col]]),
    hpo_id = str_trim(as.character(.data[[g2p_hpo_id_col]])),
    hpo_name = str_trim(as.character(.data[[g2p_hpo_name_col]])),
    frequency = str_trim(as.character(.data[[g2p_freq_col]])),
    disease_id = str_trim(as.character(.data[[g2p_disease_col]])),
    disease_mim = gsub("[^0-9]", "", disease_id)
  ) %>%
  filter(gene_id != "", hpo_id != "")

omimgenehpo <- gene2hpo %>%
  filter(gene_id %in% target_genes | disease_mim %in% target_mims)

if (nrow(omimgenehpo) == 0) {
  warning("No HPO gene-phenotype annotation matched target OMIM genes/MIM IDs.")
}

message2("HPO rows matched to target genes/diseases: ", nrow(omimgenehpo))

omimgenehpo <- omimgenehpo %>%
  mutate(
    hpo_id = ifelse(hpo_id == "HP:0025356", "HP:0001263", hpo_id),
    level2hpo = vapply(hpo_id, function(term) paste(get_level2_ancestors(term, hpo, level2_terms), collapse = ","), character(1))
  )

gene_hpo_summary <- omimgenehpo %>%
  group_by(gene_id) %>%
  summarise(
    gene_hpo_ids = paste(sort(unique(hpo_id)), collapse = ","),
    gene_hpo_ids2 = paste(sort(unique(c(split_ids(paste(hpo_id, collapse = ",")), split_ids(paste(level2hpo, collapse = ","))))), collapse = ","),
    .groups = "drop"
  )

omimgenehpo2 <- omim_core %>%
  left_join(gene_hpo_summary, by = "gene_id") %>%
  mutate(
    gene_hpo_ids = ifelse(is.na(gene_hpo_ids), "", gene_hpo_ids),
    gene_hpo_ids2 = ifelse(is.na(gene_hpo_ids2), "", gene_hpo_ids2)
  )

write_tsv(omimgenehpo2, paste0(out_prefix, ".target_gene_hpo_summary.tsv"))

# -----------------------------
# SNV table
# -----------------------------
message2("Reading SNV table: ", snv_file)
snv <- safe_read_tsv(snv_file)
input_snv_rows <- nrow(snv)

snv <- rename_if_present(snv, "Sample", "SampleID")
snv <- rename_if_present(snv, "Matched_OMIM_gene", "gene_id")
snv <- rename_if_present(snv, "Inhouse_negative_freq", "freq.in.neg4128")
snv <- rename_if_present(snv, "Inhouse_positive_freq", "freq.in.pos2035")

if (!("SampleID" %in% colnames(snv)) && sample_id_col %in% colnames(snv)) {
  colnames(snv)[colnames(snv) == sample_id_col] <- "SampleID"
}
if (!("gene_id" %in% colnames(snv))) {
  if ("Gene_refGene" %in% colnames(snv)) {
    snv$gene_id <- snv$Gene_refGene
  } else {
    stop("SNV table must contain gene_id, Matched_OMIM_gene, or Gene_refGene", call. = FALSE)
  }
}

required_snv_cols <- c("SampleID", "Chr", "Pos", "Ref", "Alt", "gene_id", "Func_refGene", "ExonicFunc_refGene")
missing_snv_cols <- setdiff(required_snv_cols, colnames(snv))
if (length(missing_snv_cols) > 0) {
  stop("Missing required SNV columns: ", paste(missing_snv_cols, collapse = ", "), call. = FALSE)
}

choose_freq_col <- function(df, opt, candidates, default_value = NA_real_) {
  if (!is.null(opt) && opt != "auto") {
    if (!(opt %in% colnames(df))) stop("Frequency column not found: ", opt, call. = FALSE)
    return(opt)
  }
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[[1]])
  NULL
}

negative_freq_col <- choose_freq_col(snv, negative_freq_col_opt, c("freq.in.neg4128", "Inhouse_negative_freq", "Inhouse_AF", "Inhouse_carrier_freq"))
positive_freq_col <- choose_freq_col(snv, positive_freq_col_opt, c("freq.in.pos2035", "Inhouse_positive_freq"))

if (is.null(negative_freq_col)) {
  warning("No negative/in-house frequency column found. Setting negative frequency to 0 for all variants.")
  snv$freq.in.neg4128 <- 0
  negative_freq_col <- "freq.in.neg4128"
}
if (is.null(positive_freq_col)) {
  warning("No positive frequency column found. Setting positive frequency to 0 for all variants.")
  snv$freq.in.pos2035 <- 0
  positive_freq_col <- "freq.in.pos2035"
}
if (!(global_af_col %in% colnames(snv))) {
  warning("Global AF column not found: ", global_af_col, ". Setting AF to NA.")
  snv[[global_af_col]] <- NA_character_
}

if (!("GT" %in% colnames(snv))) {
  snv$GT <- if ("Info" %in% colnames(snv)) vapply(snv$Info, parse_gt, character(1)) else "."
}
if (!("Info" %in% colnames(snv))) snv$Info <- "."
if (!("pLi_refGene" %in% colnames(snv))) snv$pLi_refGene <- NA
if (!("REVEL_score" %in% colnames(snv))) snv$REVEL_score <- NA
if (!("CADD16_phread_dbnsfp" %in% colnames(snv))) snv$CADD16_phread_dbnsfp <- NA
if (!("hgvsc" %in% colnames(snv))) snv$hgvsc <- if ("AAChange_refGene" %in% colnames(snv)) snv$AAChange_refGene else "."
if (!("LOF" %in% colnames(snv))) snv$LOF <- "."

snv <- snv %>%
  mutate(
    SampleID = clean_id(SampleID),
    Chr = gsub("^chr", "", as.character(Chr), ignore.case = TRUE),
    Pos = as.character(Pos),
    Ref = as.character(Ref),
    Alt = as.character(Alt),
    gene_id = clean_id(gene_id),
    snv = paste(Chr, Pos, Ref, Alt, sep = ":"),
    AF_num = to_num(.data[[global_af_col]]),
    freq.in.neg4128 = to_num(.data[[negative_freq_col]]),
    freq.in.pos2035 = to_num(.data[[positive_freq_col]]),
    pLi_num = to_num(pLi_refGene),
    REVEL_num = to_num(REVEL_score),
    CADD_num = to_num(CADD16_phread_dbnsfp),
    is_lof = is_lof_variant(Func_refGene, ExonicFunc_refGene, LOF),
    is_protein_altering = is_protein_altering_variant(Func_refGene, ExonicFunc_refGene)
  ) %>%
  distinct(SampleID, snv, gene_id, .keep_all = TRUE)

noncoding_exclude <- c("intergenic", "intronic", "downstream", "upstream;downstream", "upstream", "UTR3", "UTR5", "UTR5;UTR3")
snv <- snv %>% filter(!(Func_refGene %in% noncoding_exclude))

if (!is.null(negative_id_file) && file.exists(negative_id_file)) {
  negative_ids <- clean_id(readLines(negative_id_file, warn = FALSE))
  negative_ids <- negative_ids[negative_ids != ""]
  snv <- snv %>% filter(SampleID %in% negative_ids)
  message2("Restricted to sample IDs in ", negative_id_file, ": ", length(unique(snv$SampleID)), " samples with candidate rows")
} else if (!is.null(negative_id_file)) {
  warning("negative-id file not found; no sample restriction applied: ", negative_id_file)
}

snv <- snv %>%
  select(-any_of(c("OMIM_Phenotype", "OMIM_Inheritance", "OMIM_MIM_number", "OMIM_Update_type", "OMIM_update_time")))

snv2 <- snv %>% inner_join(omimgenehpo2, by = "gene_id")
message2("SNV rows after target gene merge: ", nrow(snv2))

snv3 <- snv2 %>%
  filter(
    is.na(freq.in.neg4128) | freq.in.neg4128 < ar_freq_cutoff,
    is.na(freq.in.pos2035) | freq.in.pos2035 < ar_freq_cutoff,
    is.na(AF_num) | AF_num < global_af_cutoff
  )
message2("SNV rows after AR/global frequency prefilter: ", nrow(snv3))

snv4 <- snv3 %>% filter(is_lof | is_protein_altering)
message2("Protein-altering/LOF/splicing rows: ", nrow(snv4))

# -----------------------------
# Patient sample HPO information
# -----------------------------
message2("Reading patient sampleinfo: ", sampleinfo_file)
sampleinfo <- safe_read_tsv(sampleinfo_file)
sampleinfo <- rename_if_present(sampleinfo, "Sample", "SampleID")
if (!("SampleID" %in% colnames(sampleinfo))) {
  stop("sampleinfo must contain SampleID or Sample column", call. = FALSE)
}
sampleinfo$SampleID <- clean_id(sampleinfo$SampleID)

if ("sample_hpo_ids" %in% colnames(sampleinfo)) {
  meta_cols <- intersect(c("SampleID", "Sex", "Age.at.diagnosis", "Age.at.diagnosis.months.", "Diagnosis.Type", "HPOTerms_Count"), colnames(sampleinfo))
  sampleinfo2 <- sampleinfo %>%
    group_by(SampleID) %>%
    summarise(
      across(all_of(setdiff(meta_cols, "SampleID")), ~ {y <- na.omit(.x); if (length(y) == 0) NA else y[[1]]}),
      sample_hpo_ids = paste(sort(unique(unlist(lapply(sample_hpo_ids, split_ids)))), collapse = ","),
      .groups = "drop"
    )
} else {
  hpo_col <- first_existing(sampleinfo, c("HPOTerms", "HPO", "HPO_ID", "hpo_id"), required = TRUE, what = "sample HPO column")
  level2_col <- first_existing(sampleinfo, c("level2hpoid", "level2_hpo_id", "level2hpo"), required = FALSE)
  meta_cols <- intersect(c("SampleID", "Sex", "Age.at.diagnosis", "Age.at.diagnosis.months.", "Diagnosis.Type", "HPOTerms_Count"), colnames(sampleinfo))
  sampleinfo2 <- sampleinfo %>%
    mutate(
      hpo_tmp = as.character(.data[[hpo_col]]),
      level2_tmp = if (!is.null(level2_col)) as.character(.data[[level2_col]]) else ""
    ) %>%
    group_by(SampleID) %>%
    summarise(
      across(all_of(setdiff(meta_cols, "SampleID")), ~ {y <- na.omit(.x); if (length(y) == 0) NA else y[[1]]}),
      sample_hpo_ids = paste(sort(unique(unlist(lapply(paste(hpo_tmp, level2_tmp, sep = ","), split_ids)))), collapse = ","),
      .groups = "drop"
    )
}

snv5 <- snv4 %>% inner_join(sampleinfo2, by = "SampleID")
message2("Rows after merging patient HPO annotations: ", nrow(snv5))

snv5 <- snv5 %>%
  mutate(
    intersecthpocount = mapply(intersect_count, sample_hpo_ids, gene_hpo_ids),
    intersecthpolist = mapply(intersect_list, sample_hpo_ids, gene_hpo_ids),
    intersecthpocount2 = mapply(intersect_count, sample_hpo_ids, gene_hpo_ids2),
    intersecthpolist2 = mapply(intersect_list, sample_hpo_ids, gene_hpo_ids2)
  )

if (match_level == "exact") {
  snv6 <- snv5 %>% filter(intersecthpocount >= min_hpo_intersect)
} else {
  snv6 <- snv5 %>% filter(intersecthpocount2 >= min_hpo_intersect)
}

snv6 <- snv6 %>%
  group_by(SampleID, snv, gene_id) %>%
  arrange(desc(intersecthpocount), desc(intersecthpocount2), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

message2("Phenotype-matched rows retained: ", nrow(snv6))
write_tsv(snv6, paste0(out_prefix, ".phenotype_matched_variants.tsv"))

# -----------------------------
# Candidate filtering by inheritance model
# -----------------------------
snv6 <- snv6 %>%
  mutate(
    is_AD_XLD = vapply(OMIM_Inheritance, has_inheritance, logical(1), codes = c("AD", "XLD")),
    is_AR_XLR = vapply(OMIM_Inheritance, has_inheritance, logical(1), codes = c("AR", "XLR")),
    hom_or_hemi_alt = vapply(GT, is_hom_or_hemi_alt, logical(1)),
    freq.in.neg4128 = ifelse(is.na(freq.in.neg4128), 0, freq.in.neg4128),
    freq.in.pos2035 = ifelse(is.na(freq.in.pos2035), 0, freq.in.pos2035)
  )

ad_neg_cutoff <- max(ad_freq_cutoff, 2 / negative_denominator)
ADbase <- snv6 %>%
  filter(
    is_AD_XLD,
    freq.in.neg4128 <= ad_neg_cutoff,
    freq.in.pos2035 == 0,
    is.na(AF_num) | AF_num < ad_freq_cutoff
  )

ADLOF <- ADbase %>% filter(!is.na(pLi_num), pLi_num >= pli_cutoff, is_lof)
ADmissense <- ADbase %>%
  filter(
    is_protein_altering,
    (
      (!is.na(REVEL_num) & !is.na(CADD_num) & REVEL_num > revel_cutoff & CADD_num > cadd_cutoff) |
        Func_refGene %in% c("ncRNA_exonic", "splicing")
    )
  )
ADsnv2 <- bind_rows(ADLOF, ADmissense) %>% distinct(SampleID, snv, gene_id, .keep_all = TRUE)

ARbase <- snv6 %>% filter(is_AR_XLR)
ARLOF <- ARbase %>% filter(!is.na(pLi_num), pLi_num >= pli_cutoff, is_lof)
ARmissense <- ARbase %>%
  filter(
    is_protein_altering,
    (is.na(REVEL_num) | REVEL_num > revel_cutoff),
    (is.na(CADD_num) | CADD_num > cadd_cutoff)
  )
ARsnv2 <- bind_rows(ARLOF, ARmissense) %>% distinct(SampleID, snv, gene_id, .keep_all = TRUE)
ARselect <- ARsnv2 %>%
  group_by(SampleID, gene_id) %>%
  mutate(gene.snv.count = n_distinct(snv)) %>%
  ungroup() %>%
  filter(gene.snv.count > 1 | hom_or_hemi_alt) %>%
  distinct(SampleID, snv, gene_id, .keep_all = TRUE)

otherBase <- snv6 %>% filter(!is_AD_XLD & !is_AR_XLR)
otherLOF <- otherBase %>% filter(!is.na(pLi_num), pLi_num >= pli_cutoff, is_lof)
otherMissense <- otherBase %>%
  filter(
    is_protein_altering,
    (is.na(REVEL_num) | REVEL_num > revel_cutoff),
    (is.na(CADD_num) | CADD_num > cadd_cutoff)
  )
otherIHsnv2 <- bind_rows(otherLOF, otherMissense) %>%
  distinct(SampleID, snv, gene_id, .keep_all = TRUE) %>%
  group_by(SampleID, Gene_refGene, gene_id) %>%
  mutate(gene.snv.count = n_distinct(snv)) %>%
  ungroup()

# Add sample HPO term names
ADsnv2 <- ADsnv2 %>% mutate(sample_hpo_terms = vapply(sample_hpo_ids, hpo_ids_to_terms, character(1), hpo_name_map = hpo_name_map))
ARselect <- ARselect %>% mutate(sample_hpo_terms = vapply(sample_hpo_ids, hpo_ids_to_terms, character(1), hpo_name_map = hpo_name_map))
otherIHsnv2 <- otherIHsnv2 %>% mutate(sample_hpo_terms = vapply(sample_hpo_ids, hpo_ids_to_terms, character(1), hpo_name_map = hpo_name_map))

# Optional clinical registry merge
merge_clinical <- function(df, pheno_xlsx, clinical_sample_col) {
  if (is.null(pheno_xlsx) || pheno_xlsx == "" || !file.exists(pheno_xlsx)) return(df)
  pheno_raw <- readxl::read_excel(pheno_xlsx)
  pheno_raw <- as.data.frame(pheno_raw, stringsAsFactors = FALSE)
  if (!(clinical_sample_col %in% colnames(pheno_raw))) {
    warning("Clinical sample column not found in Excel: ", clinical_sample_col, ". Skip clinical merge.")
    return(df)
  }
  colnames(pheno_raw)[colnames(pheno_raw) == clinical_sample_col] <- "SampleID"
  pheno_raw$SampleID <- clean_id(pheno_raw$SampleID)
  keep <- intersect(c("SampleID"), colnames(pheno_raw))
  pheno_keep <- pheno_raw %>% select(all_of(keep)) %>% distinct(SampleID, .keep_all = TRUE)
  left_join(df, pheno_keep, by = "SampleID", suffix = c("", ".clinical"))
}

ADsnv2_pheno <- merge_clinical(ADsnv2, pheno_xlsx, clinical_sample_col)
ARselect_pheno <- merge_clinical(ARselect, pheno_xlsx, clinical_sample_col)
otherIHsnv2_pheno <- merge_clinical(otherIHsnv2, pheno_xlsx, clinical_sample_col)

# -----------------------------
# Outputs
# -----------------------------
write_tsv(ADsnv2_pheno, paste0(out_prefix, ".AD_XLD_candidates.tsv"))
write_tsv(ARselect_pheno, paste0(out_prefix, ".AR_XLR_candidates.tsv"))
write_tsv(otherIHsnv2_pheno, paste0(out_prefix, ".other_inheritance_candidates.tsv"))

write_excel_if_possible(
  list(
    AD_XLD = ADsnv2_pheno,
    AR_XLR = ARselect_pheno,
    other_inheritance = otherIHsnv2_pheno,
    phenotype_matched = snv6
  ),
  paste0(out_prefix, ".candidate_variants.xlsx")
)

summary_tbl <- data.frame(
  metric = c(
    "input_snv_rows",
    "target_omim_genes",
    "target_gene_hpo_rows",
    "snv_rows_after_target_gene_merge",
    "snv_rows_after_frequency_prefilter",
    "protein_altering_or_lof_rows",
    "rows_after_sample_hpo_merge",
    "phenotype_matched_rows",
    "phenotype_matched_unique_sample_variant",
    "AD_XLD_candidate_variants",
    "AD_XLD_candidate_samples",
    "AD_XLD_candidate_genes",
    "AR_XLR_candidate_variants",
    "AR_XLR_candidate_samples",
    "AR_XLR_candidate_genes",
    "other_candidate_variants",
    "other_candidate_samples",
    "other_candidate_genes"
  ),
  value = c(
    input_snv_rows,
    length(target_genes),
    nrow(omimgenehpo),
    nrow(snv2),
    nrow(snv3),
    nrow(snv4),
    nrow(snv5),
    nrow(snv6),
    nrow(distinct(snv6, SampleID, snv)),
    nrow(distinct(ADsnv2_pheno, SampleID, snv)),
    length(unique(ADsnv2_pheno$SampleID)),
    length(unique(ADsnv2_pheno$gene_id)),
    nrow(distinct(ARselect_pheno, SampleID, snv)),
    length(unique(ARselect_pheno$SampleID)),
    length(unique(ARselect_pheno$gene_id)),
    nrow(distinct(otherIHsnv2_pheno, SampleID, snv)),
    length(unique(otherIHsnv2_pheno$SampleID)),
    length(unique(otherIHsnv2_pheno$gene_id))
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary_tbl, paste0(out_prefix, ".qc_summary.tsv"))
write_excel_if_possible(list(summary = summary_tbl), paste0(out_prefix, ".qc_summary.xlsx"))

message2("Done.")
message2("AD/XLD candidates: ", nrow(ADsnv2_pheno), " rows; ", length(unique(ADsnv2_pheno$SampleID)), " samples")
message2("AR/XLR candidates: ", nrow(ARselect_pheno), " rows; ", length(unique(ARselect_pheno$SampleID)), " samples")
message2("Other inheritance candidates: ", nrow(otherIHsnv2_pheno), " rows; ", length(unique(otherIHsnv2_pheno$SampleID)), " samples")
message2("Output prefix: ", out_prefix)
