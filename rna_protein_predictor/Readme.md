# Quantitative RNA-to-protein predictability across cancers

This repository contains a prespecified, cross-validated analysis of how well
matched RNA abundance predicts quantitative protein abundance across human
cancers.

The primary scientific question is:

> Which gene properties determine how faithfully RNA abundance predicts protein
> abundance across cancers?

The active project uses matched CPTAC transcriptomics and mass-spectrometry
proteomics from BRCA, COAD, GBM, and PDAC. Earlier HPA/Xena classification and
assay-disagreement experiments were removed from the active tree after the
quantitative project was frozen. They remain recoverable from the Git tag
`pre-cptac-streamline-20260828`.

## Current status

Completed checkpoints:

1. Download and match four CPTAC cohorts.
2. Evaluate RNA-based protein prediction in held-out patients.
3. Test OLS, Ridge, and Huber estimator sensitivity.
4. Calculate patient- and gene-bootstrap confidence intervals.
5. Evaluate cross-cancer stability and held-out-cancer transfer.
6. Build an 8,681-gene cross-cancer predictability atlas.
7. Fit a repeated held-out-gene technical baseline.

The quantitative checkpoint is **PASS**. RNA improves held-out protein
prediction in all four cancers, and basic measurement characteristics explain
about one-third of cross-gene predictability variation.

The frozen design and go/no-go rules are documented in
[`docs/QUANTITATIVE_CPTAC_BENCHMARK.md`](docs/QUANTITATIVE_CPTAC_BENCHMARK.md).

## Validated findings

### Held-out-patient prediction

WashU RNA FPKM values are transformed as `log2(FPKM + 1)`. Within each cancer,
five-fold patient-level cross-validation compares a training-fold gene-mean
protein predictor with a gene-specific RNA model. Patient folds are shared
across genes, and each gene requires at least 30 paired observations.

| Cancer | Tumors | Eligible genes | Gene-mean RMSE | RNA RMSE | Relative MSE improvement |
|---|---:|---:|---:|---:|---:|
| BRCA | 119 | 10,256 | 0.5597 | 0.4693 | 29.7% |
| COAD | 96 | 7,261 | 0.4445 | 0.4031 | 17.8% |
| GBM | 99 | 10,443 | 0.4877 | 0.4045 | 31.2% |
| PDAC | 105 | 9,315 | 0.4719 | 0.4193 | 21.0% |

RNA improves gene-level MSE for 74.5%-87.1% of eligible genes. Excluding the ten
most harmful extrapolating genes in a sensitivity audit changes aggregate
improvement by only 0.3-2.5 percentage points.

### Estimator and bootstrap robustness

OLS, nested RidgeCV, and Huber regression use identical outer patient folds and
produce similar cohort-level improvements.

| Cancer | OLS improvement | Patient-bootstrap 95% CI | Ridge | Huber |
|---|---:|---:|---:|---:|
| BRCA | 29.7% | 26.3%-32.7% | 29.7% | 29.5% |
| COAD | 17.8% | 15.0%-20.4% | 17.8% | 17.7% |
| GBM | 31.2% | 27.5%-34.4% | 31.4% | 30.6% |
| PDAC | 21.0% | 18.4%-23.8% | 21.3% | 21.0% |

Gene-bootstrap intervals also show positive median improvement in every cancer.
Fewer than 0.3% of predictions fall outside the corresponding training protein
range. Isolated extreme extrapolations are retained and documented.

### Cross-cancer stability

- Pairwise gene-profile Spearman correlations: 0.504-0.634.
- Median held-out-cancer correlations: 0.318-0.472.
- Positive held-out-cancer direction: 92.4%-96.3% of evaluated genes.

The held-out-cancer analysis evaluates transfer of coupling direction and
patient ranking. It does not establish absolute calibration on a new
proteomics intensity scale.

### Gene-predictability atlas

The primary gene-level outcome is median cross-cancer Ridge MSE improvement
relative to the gene-mean baseline. Genes must be evaluated in at least three
cancers.

| Atlas property | Result |
|---|---:|
| Eligible genes | 8,681 |
| Observed in three cancers | 2,007 |
| Observed in four cancers | 6,674 |
| Median improvement | 14.4% |
| Positive median improvement | 88.4% |
| Positive in every observed cancer | 55.2% |

The primary explanatory analysis retains this continuous outcome rather than
creating an arbitrary high/low gene split.

### Held-out-gene technical baseline

A repeated five-fold Ridge model predicts gene-level performance from RNA and
protein abundance, variability, protein coverage, paired-patient count, and
cancer coverage.

| Metric | Technical model | 95% gene-bootstrap CI |
|---|---:|---:|
| Cross-validated R-squared | 0.330 | 0.309-0.351 |
| Spearman correlation | 0.607 | 0.592-0.621 |
| RMSE | 0.146 | 0.144-0.149 |
| MAE | 0.115 | 0.113-0.117 |

RNA variability is the strongest technical correlate. Individual standardized
coefficients are not treated as causal because several abundance and
variability features are correlated.

## Data source

[Clinical Proteomic Tumor Analysis Consortium (CPTAC)](https://proteomic.datacommons.cancer.gov/pdc/cptac-pancancer)
matched transcriptomics and mass-spectrometry proteomics are accessed through
the Python `cptac` package.

Discovery inputs:

- WashU transcriptomics.
- BCM reference-intensity-normalized proteomics.
- BRCA, COAD, GBM, and PDAC tumors.

PDAC contains 15 normal samples that are preserved locally but excluded from
the primary tumor benchmark.

## Installation

Python 3.9 or newer is required.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

## Reproduce the analysis

Run commands from `rna_protein_predictor/`.

### 1. Download and match cohorts

```bash
python3 scripts/cptac_pilot.py --cancer coad --rna-source washu --protein-source bcm
python3 scripts/cptac_pilot.py --cancer brca --rna-source washu --protein-source bcm
python3 scripts/cptac_pilot.py --cancer gbm  --rna-source washu --protein-source bcm
python3 scripts/cptac_pilot.py --cancer pdac --rna-source washu --protein-source bcm
```

### 2. Run the quantitative benchmark

```bash
python3 scripts/benchmark_cptac_quantitative.py
```

### 3. Validate models and uncertainty

```bash
python3 scripts/validate_cptac_quantitative_models.py \
  --n-jobs 4 \
  --bootstrap 2000
```

### 4. Build the gene atlas

```bash
python3 scripts/build_gene_predictability_atlas.py
```

### 5. Fit the technical baseline

```bash
python3 scripts/model_gene_predictability_technical.py \
  --folds 5 \
  --repeats 20 \
  --bootstrap 2000
```

## Active scripts

| Script | Purpose |
|---|---|
| `scripts/cptac_pilot.py` | Download, match, and audit one CPTAC cohort |
| `scripts/benchmark_cptac_quantitative.py` | Primary held-out-patient benchmark and cross-cancer stability |
| `scripts/validate_cptac_quantitative_models.py` | OLS/Ridge/Huber validation and bootstrap intervals |
| `scripts/build_gene_predictability_atlas.py` | Continuous cross-cancer gene-predictability atlas |
| `scripts/model_gene_predictability_technical.py` | Repeated held-out-gene technical baseline |

## Important outputs

| Output | Purpose |
|---|---|
| `outputs/reports/cptac_quantitative_benchmark.json` | Primary benchmark summary |
| `outputs/reports/cptac_quantitative_model_validation.json` | Robust-model validation |
| `outputs/reports/quantitative_cptac_checkpoint.json` | Initial frozen go/no-go decision |
| `outputs/reports/quantitative_model_robustness_checkpoint.json` | Robustness checkpoint |
| `outputs/reports/cptac_gene_predictability_atlas.json` | Atlas audit |
| `outputs/reports/gene_predictability_technical_model.json` | Technical baseline |
| `outputs/tables/cptac_quantitative_cv_summary.csv` | Cancer benchmark metrics |
| `outputs/tables/cptac_quantitative_bootstrap_intervals.csv` | Bootstrap intervals |
| `outputs/tables/cptac_coupling_pairwise_stability.csv` | Cross-cancer stability |
| `outputs/tables/cptac_quantitative_loco_transfer.csv` | Held-out-cancer transfer |
| `outputs/tables/cptac_gene_predictability_atlas.csv` | Primary gene-level atlas |
| `outputs/tables/gene_predictability_technical_bootstrap_intervals.csv` | Technical-model uncertainty |

Downloaded matrices and generated tables are ignored by Git. Compact reports
and checkpoint summaries may be versioned.

## Repository structure

```text
rna_protein_predictor/
├── data/
│   └── cptac/              # Local matched cohort matrices; ignored by Git
├── docs/
│   └── QUANTITATIVE_CPTAC_BENCHMARK.md
├── outputs/
│   ├── reports/            # Compact result summaries
│   ├── state/              # Local checkpoints; ignored by Git
│   └── tables/             # Generated tables; ignored by Git
├── scripts/                # Five active pipeline entry points
├── .gitignore
├── Readme.md
└── requirements.txt
```

## Interpretation boundaries

- The models predict relative protein abundance within cohorts, not absolute
  calibration across proteomics studies.
- RNA-protein association does not establish direct translational regulation.
- Bulk tumor purity and cell composition can influence both assays.
- Protein missingness and dynamic range influence apparent predictability and
  are included in the technical baseline.
- Alternative processing workflows on the same patients provide technical
  robustness, not independent biological replication.
- Rare out-of-range predictions are retained rather than removed post hoc.

## Remaining publication checkpoints

1. Retrieve and freeze independent UniProt/GO biological annotations.
2. Require at least 70% annotation coverage of the 8,681-gene atlas.
3. Compare combined technical/biological models with the frozen technical
   baseline using identical held-out-gene folds.
4. Require a paired-bootstrap interval excluding zero, permutation
   `P < 0.05`, and incremental cross-validated R-squared of at least 0.02 for
   a meaningful biological contribution.
5. Validate with an alternative proteomics workflow and, where possible, an
   independent biological cohort.
6. Evaluate tumor purity and composition sensitivity without post hoc sample
   removal.
7. Generate final figures only after biological and external validation.

## Recovery of legacy work

The legacy HPA/Xena, binary-classification, discordance, localization, and
cross-assay analyses were removed from the active tree during the CPTAC
streamlining checkpoint. Their committed versions remain available from:

```bash
git show pre-cptac-streamline-20260828
```

To inspect the old project without changing the current branch:

```bash
git worktree add ../rna_protein_predictor_legacy \
  pre-cptac-streamline-20260828
```

Local files moved during cleanup remain recoverable from:

```text
/Users/kj525/.Trash/rna_protein_predictor_legacy_20260828/
```

## Citation

Data were generated by the Clinical Proteomic Tumor Analysis Consortium
(NCI/NIH). A manuscript must cite the exact CPTAC cohort datasets, processing
sources, package version, and cohort-specific publications used by the frozen
analysis.
