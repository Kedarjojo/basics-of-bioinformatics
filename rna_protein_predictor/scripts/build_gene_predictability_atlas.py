#!/usr/bin/env python3
"""Build the gene-level CPTAC predictability atlas for checkpoint 4."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    from benchmark_cptac_quantitative import DEFAULT_CANCERS, load_cancer
except ModuleNotFoundError:  # Supports module-style execution in development.
    from src.benchmark_cptac_quantitative import DEFAULT_CANCERS, load_cancer


ROOT = Path(__file__).resolve().parents[1]
PREDICTIONS = ROOT / "outputs" / "tables" / "cptac_quantitative_model_predictions.csv.gz"


def _safe_spearman(x: pd.Series, y: pd.Series) -> float:
    valid = x.notna() & y.notna()
    if valid.sum() < 3 or x[valid].nunique() < 2 or y[valid].nunique() < 2:
        return np.nan
    return float(spearmanr(x[valid], y[valid]).statistic)


def prediction_metrics(pred: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for (cancer, gene), group in pred.groupby(["cancer", "gene"], sort=True):
        baseline_sse = float(group["se_gene_mean"].sum())
        row: dict[str, object] = {
            "cancer": cancer,
            "gene": gene,
            "n_oof_pairs": len(group),
            "baseline_sse": baseline_sse,
        }
        for model in ("ols", "ridge", "huber"):
            model_sse = float(group[f"se_{model}"].sum())
            row[f"{model}_sse"] = model_sse
            row[f"{model}_mse_improvement"] = (
                1 - model_sse / baseline_sse if baseline_sse > 0 else np.nan
            )
            row[f"{model}_prediction_spearman"] = _safe_spearman(
                group["observed"], group[f"pred_{model}"]
            )
            row[f"{model}_outside_range_fraction"] = float(
                group[f"outside_train_{model}"].mean()
            )
        rows.append(row)
    return pd.DataFrame(rows)


def technical_metrics(
    cancers: list[str], min_pairs: int
) -> pd.DataFrame:
    rows: list[dict] = []
    for cancer in cancers:
        rna, protein = load_cancer(cancer)
        for gene in rna.columns:
            x = rna[gene]
            y = protein[gene]
            paired = x.notna() & y.notna()
            if paired.sum() < min_pairs:
                continue
            rows.append(
                {
                    "cancer": cancer.upper(),
                    "gene": gene,
                    "matched_samples": len(rna),
                    "n_paired": int(paired.sum()),
                    "protein_coverage": float(y.notna().mean()),
                    "rna_missing_fraction": float(x.isna().mean()),
                    "rna_mean": float(x[paired].mean()),
                    "rna_sd": float(x[paired].std(ddof=1)),
                    "rna_iqr": float(x[paired].quantile(0.75) - x[paired].quantile(0.25)),
                    "protein_mean": float(y[paired].mean()),
                    "protein_sd": float(y[paired].std(ddof=1)),
                    "protein_iqr": float(y[paired].quantile(0.75) - y[paired].quantile(0.25)),
                    "observed_rna_protein_spearman": _safe_spearman(x, y),
                }
            )
    return pd.DataFrame(rows)


def _iqr(series: pd.Series) -> float:
    return float(series.quantile(0.75) - series.quantile(0.25))


def aggregate_gene_atlas(per_cancer: pd.DataFrame, min_cancers: int) -> pd.DataFrame:
    rows: list[dict] = []
    for gene, group in per_cancer.groupby("gene", sort=True):
        if group["cancer"].nunique() < min_cancers:
            continue
        ridge = group["ridge_mse_improvement"]
        rows.append(
            {
                "gene": gene,
                "cancers_observed": int(group["cancer"].nunique()),
                "cancers": ";".join(sorted(group["cancer"].unique())),
                "primary_median_ridge_improvement": float(ridge.median()),
                "ridge_improvement_min": float(ridge.min()),
                "ridge_improvement_max": float(ridge.max()),
                "ridge_improvement_iqr": _iqr(ridge),
                "positive_ridge_cancers": int((ridge > 0).sum()),
                "median_ols_improvement": float(group["ols_mse_improvement"].median()),
                "median_huber_improvement": float(group["huber_mse_improvement"].median()),
                "median_ridge_prediction_spearman": float(
                    group["ridge_prediction_spearman"].median()
                ),
                "median_observed_rna_protein_spearman": float(
                    group["observed_rna_protein_spearman"].median()
                ),
                "median_rna_mean": float(group["rna_mean"].median()),
                "median_rna_sd": float(group["rna_sd"].median()),
                "median_rna_iqr": float(group["rna_iqr"].median()),
                "median_protein_mean": float(group["protein_mean"].median()),
                "median_protein_sd": float(group["protein_sd"].median()),
                "median_protein_iqr": float(group["protein_iqr"].median()),
                "median_protein_coverage": float(group["protein_coverage"].median()),
                "minimum_protein_coverage": float(group["protein_coverage"].min()),
                "median_paired_patients": float(group["n_paired"].median()),
                "median_ridge_outside_range_fraction": float(
                    group["ridge_outside_range_fraction"].median()
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        "primary_median_ridge_improvement", ascending=False
    )


def overview(
    pred: pd.DataFrame,
    per_cancer: pd.DataFrame,
    atlas: pd.DataFrame,
    min_cancers: int,
) -> dict:
    target = atlas["primary_median_ridge_improvement"]
    return {
        "analysis": "CPTAC gene predictability atlas",
        "primary_outcome": "median cross-cancer Ridge MSE improvement",
        "minimum_cancers": min_cancers,
        "prediction_rows": len(pred),
        "gene_cancer_rows": len(per_cancer),
        "eligible_genes": len(atlas),
        "genes_by_cancers_observed": {
            str(k): int(v)
            for k, v in atlas["cancers_observed"].value_counts().sort_index().items()
        },
        "primary_outcome_quantiles": {
            str(q): float(value)
            for q, value in target.quantile(
                [0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1]
            ).items()
        },
        "positive_primary_outcome_fraction": float((target > 0).mean()),
        "positive_in_all_observed_cancers_fraction": float(
            (atlas["positive_ridge_cancers"] == atlas["cancers_observed"]).mean()
        ),
        "complete_technical_covariates": bool(
            atlas[
                [
                    "median_rna_mean",
                    "median_rna_sd",
                    "median_protein_mean",
                    "median_protein_sd",
                    "median_protein_coverage",
                    "median_paired_patients",
                ]
            ].notna().all().all()
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancers", nargs="+", default=list(DEFAULT_CANCERS))
    parser.add_argument("--min-pairs", type=int, default=30)
    parser.add_argument("--min-cancers", type=int, default=3)
    parser.add_argument("--predictions", type=Path, default=PREDICTIONS)
    args = parser.parse_args()

    if not args.predictions.exists():
        raise FileNotFoundError(
            f"Missing model predictions: {args.predictions}. Run "
            "scripts/validate_cptac_quantitative_models.py first."
        )
    pred = pd.read_csv(args.predictions)
    required = {
        "cancer", "gene", "observed", "se_gene_mean",
        "se_ols", "se_ridge", "se_huber",
        "pred_ols", "pred_ridge", "pred_huber",
        "outside_train_ols", "outside_train_ridge", "outside_train_huber",
    }
    missing = sorted(required.difference(pred.columns))
    if missing:
        raise ValueError(f"Prediction table is missing columns: {missing}")

    performance = prediction_metrics(pred)
    technical = technical_metrics(args.cancers, args.min_pairs)
    per_cancer = performance.merge(
        technical,
        on=["cancer", "gene"],
        how="inner",
        validate="one_to_one",
    )
    atlas = aggregate_gene_atlas(per_cancer, args.min_cancers)
    report = overview(pred, per_cancer, atlas, args.min_cancers)

    tables = ROOT / "outputs" / "tables"
    reports = ROOT / "outputs" / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    per_cancer.to_csv(
        tables / "cptac_gene_predictability_by_cancer.csv", index=False
    )
    atlas.to_csv(tables / "cptac_gene_predictability_atlas.csv", index=False)
    (reports / "cptac_gene_predictability_atlas.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
