# Prespecified quantitative CPTAC benchmark

## Status

This document freezes the go/no-go analysis before inspecting its results. The
binary HPA endpoint and the HPA--CPTAC disagreement atlas are not primary
outcomes in this benchmark.

## Primary question

How accurately does matched RNA abundance predict quantitative protein
abundance in unseen patients, and is gene-level RNA--protein coupling stable
across cancers?

## Cohorts

- BRCA, COAD, GBM and PAAD/PDAC.
- Tumours only for the primary analysis.
- PDAC normals are retained but excluded from model fitting.
- WashU RNA and BCM proteomics are the discovery inputs.
- A second proteomics workflow is confirmatory and must not be called an
  independent biological cohort when it contains the same patients.

## Primary estimands

1. **Held-out-patient prediction:** five-fold patient-level cross-validation
   within each cancer. Splits are shared across genes.
2. **Increment over baseline:** reduction in squared error for a gene-wise
   linear RNA model relative to a training-fold gene-mean protein predictor.
3. **Within-gene ranking:** out-of-fold Spearman correlation between predicted
   and observed protein across patients.
4. **Held-out-cancer transfer:** leave one cancer out, fit gene-specific slopes
   in standardized training-cancer data, and test whether those slopes rank
   patients correctly in the unseen cancer.
5. **Gene stability:** concordance of observed gene-level RNA--protein
   correlations across cancers, including direction consistency.

## Models in the frozen first benchmark

- `gene_mean`: training-fold mean protein abundance for each gene.
- `rna_linear`: one gene-specific ordinary least-squares model using the
  matched RNA abundance for that gene.

More complex models, CNV, purity, composition and biological annotations are
secondary analyses. They may be added only after the two frozen models have
been evaluated and must use identical patient folds.

## Inclusion rules

- At least 30 paired non-missing RNA/protein observations per gene and cancer.
- Protein values must vary in the relevant training and test data.
- No imputation of missing protein outcomes.
- RNA missingness is reported; unexpected RNA missingness is not silently
  imputed in the first benchmark.
- Genes are identified by the symbols in the matched CPTAC matrices; ambiguous
  source identifiers have already been collapsed by the pilot loader.

## Metrics

Primary:

- Pooled out-of-fold RMSE by cancer and model.
- Relative MSE improvement of `rna_linear` over `gene_mean`.
- Median gene-level out-of-fold Spearman correlation.

Secondary:

- Pearson correlation, MAE and conventional out-of-fold R-squared.
- Distribution of gene-level correlations rather than only their mean.
- Leave-one-cancer-out gene-level Spearman correlation.

The leave-one-cancer-out analysis evaluates transfer of RNA--protein coupling,
not absolute cross-study protein calibration. Protein reference-intensity
scales can differ between cohorts, so training cancers are standardized before
slope estimation and the held-out cancer is assessed using within-gene
association metrics.

## Uncertainty

Final confidence intervals will use a patient bootstrap for prediction metrics
and a gene bootstrap for distributions of gene-level performance. Candidate
gene lists are not interpreted before these intervals and stability analyses
are complete.

## Go/no-go rules

The project advances to biological-feature modelling only if all of the
following hold:

1. `rna_linear` improves pooled held-out-patient MSE over `gene_mean` in at
   least three of four cancers.
2. The improvement is not driven solely by a small number of genes.
3. Gene-level coupling shows positive cross-cancer concordance, with a
   reproducible subset whose direction is consistent in at least three
   cancers.
4. Held-out-cancer slopes provide positive median rank association in at least
   three of four held-out cancers.

If these conditions fail, we report that RNA--protein predictability is too
cohort-dependent for the proposed cross-cancer biological-property paper and
do not create another post hoc central hypothesis.

## Staged answers to the five checkpoints

| Checkpoint | Frozen analysis | Stage |
|---|---|---|
| Held-out-patient accuracy | Five-fold within-cancer CV | 1 |
| Held-out-cancer generalization | Standardized gene-wise LOCO slopes | 1 |
| Stability of gene performance | Cross-cancer concordance and direction | 1 |
| Biological explanation | Independent annotations; nested comparison | 2 |
| Independent/orthogonal validation | Alternative proteomics workflow and, where available, independent cohort | 3 |

