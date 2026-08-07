# Workflow overview

This pipeline converts periodic OMIM updates into candidate variant tables for manual clinical review.

1. Download HPO reference files (`hp.obo` and `genes_to_phenotype.txt`).
2. Collect OMIM update entries from monthly update pages, saved email files, or an IMAP mailbox.
3. Retain disease entries from `New Entries` and `New Clinical Synopses` by default.
4. Link OMIM disease MIM identifiers to gene symbols, inheritance terms, and HPO annotations using `genes_to_phenotype.txt`.
5. Extract rare variants in updated genes from initially unresolved cases.
6. Annotate in-house carrier and allele frequencies.
7. Prioritize candidate variants using inheritance model, allele frequency, functional consequence, damaging scores, gene constraint, and HPO phenotype matching.
8. Export inheritance-stratified candidate tables for manual review, Sanger validation, segregation analysis, and ACMG/AMP-based classification.

The pipeline does not perform automated clinical diagnosis. Final interpretation requires expert review.
