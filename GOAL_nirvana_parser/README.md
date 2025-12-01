# GOAL Nirvana Annotation Parser

### GitHub README File
A lightweight toolkit for converting GOAL DRAGEN Nirvana JSON output into a clean TSV with correct transcript selection.

## Description
Parses the `positions` block of Nirvana JSON, keeps only **PASS** variants, selects a single transcript per variant (**MANE-first**), and outputs a review-ready TSV.

## Features
- Parse Nirvana JSON → TSV
- Filter to **PASS** variants
- Split multi-variant positions
- Transcript selection:
  - **Primary:** MANE Select
  - **Fallback:** `mRNA` + `isCanonical == true`
- Outputs variant, transcript, ClinVar, population, and COSMIC fields

## Input
- Nirvana JSON file (`*.vcf.annotated.json.gz`) from GOAL DRAGEN

## Output

TSV with key fields:
- **Variant:** chr, pos, ref/alt, filters, depth, genotype
- **Transcript:** transcript name, consequence, hgvsc/hgvsp
- **ClinVar:** ID, significance, reviewStatus
- **Population:** gnomad/1000g/topmed AFs
- **COSMIC:** IDs

## Validation

- Unit tests for parser and transcript logic
- Manual spot-checks of representative variants
