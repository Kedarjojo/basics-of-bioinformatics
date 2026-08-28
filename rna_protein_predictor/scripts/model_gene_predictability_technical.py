#!/usr/bin/env python3
"""Held-out-gene technical baseline for the CPTAC predictability atlas."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import RepeatedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "outputs" / "tables" / "cptac_gene_predictability_atlas.csv"
TARGET = "primary_median_ridge_improvement"
FEATURES = [
    "median_rna_mean",
    "median_rna_sd",
    "median_rna_iqr",
    "median_protein_mean",
    "median_protein_sd",
    "median_protein_iqr",
    "median_protein_coverage",
    "minimum_protein_coverage",
    "median_paired_patients",
    "cancers_observed",
]
ALPHAS = np.logspace(-4, 4, 17)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(spearmanr(x, y).statistic)


def cross_validate(
    data: pd.DataFrame, folds: int, repeats: int, seed: int
) -> pd.DataFrame:
    x = data[FEATURES].to_numpy(float)
    y = data[TARGET].to_numpy(float)
    genes = data["gene"].astype(str).to_numpy()
    splitter = RepeatedKFold(
        n_splits=folds,
        n_repeats=repeats,
        random_state=seed,
    )
    rows: list[dict] = []
    for split_index, (train, test) in enumerate(splitter.split(x)):
        repeat = split_index // folds
        fold = split_index % folds
        intercept_only_prediction = np.full(len(test), np.mean(y[train]))
        model = make_pipeline(
            StandardScaler(),
            RidgeCV(alphas=ALPHAS),
        ).fit(x[train], y[train])
        technical_prediction = model.predict(x[test])
        alpha = float(model[-1].alpha_)
        for position, intercept_only_value, technical_value in zip(
            test, intercept_only_prediction, technical_prediction
        ):
            rows.append(
                {
                    "gene": genes[position],
                    "repeat": repeat,
                    "fold": fold,
                    "observed": y[position],
                    "pred_intercept_only": float(intercept_only_value),
                    "pred_technical": float(technical_value),
                    "ridge_alpha": alpha,
                }
            )
    return pd.DataFrame(rows)


def repeat_metrics(pred: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for repeat, group in pred.groupby("repeat", sort=True):
        observed = group["observed"].to_numpy(float)
        for model in ("intercept_only", "technical"):
            fitted = group[f"pred_{model}"].to_numpy(float)
            rows.append(
                {
                    "repeat": int(repeat),
                    "model": model,
                    "r2": float(r2_score(observed, fitted)),
                    "rmse": float(np.sqrt(mean_squared_error(observed, fitted))),
                    "mae": float(mean_absolute_error(observed, fitted)),
                    "spearman": safe_spearman(observed, fitted),
                }
            )
    return pd.DataFrame(rows)


def bootstrap_average_predictions(
    pred: pd.DataFrame, replicates: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    averaged = (
        pred.groupby("gene", sort=True)
        .agg(
            observed=("observed", "first"),
            pred_intercept_only=("pred_intercept_only", "mean"),
            pred_technical=("pred_technical", "mean"),
        )
        .reset_index()
    )
    rng = np.random.default_rng(seed)
    rows = []
    n = len(averaged)
    for replicate in range(replicates):
        draw = averaged.iloc[rng.integers(0, n, n)]
        observed = draw["observed"].to_numpy(float)
        for model in ("intercept_only", "technical"):
            fitted = draw[f"pred_{model}"].to_numpy(float)
            rows.append(
                {
                    "replicate": replicate,
                    "model": model,
                    "r2": float(r2_score(observed, fitted)),
                    "rmse": float(np.sqrt(mean_squared_error(observed, fitted))),
                    "mae": float(mean_absolute_error(observed, fitted)),
                    "spearman": safe_spearman(observed, fitted),
                }
            )
    return averaged, pd.DataFrame(rows)


def bootstrap_intervals(bootstrap: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, group in bootstrap.groupby("model", sort=True):
        for metric in ("r2", "rmse", "mae", "spearman"):
            values = group[metric].dropna().to_numpy(float)
            rows.append(
                {
                    "model": model,
                    "metric": metric,
                    "estimate": float(np.median(values)),
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                }
            )
    return pd.DataFrame(rows)


def full_model_coefficients(data: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    x = data[FEATURES].to_numpy(float)
    y = data[TARGET].to_numpy(float)
    model = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS)).fit(x, y)
    coefficients = pd.DataFrame(
        {
            "feature": FEATURES,
            "standardized_coefficient": model[-1].coef_,
        }
    ).sort_values("standardized_coefficient", ascending=False)
    return coefficients, float(model[-1].alpha_)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, default=ATLAS)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260828)
    args = parser.parse_args()

    data = pd.read_csv(args.atlas)
    required = {"gene", TARGET, *FEATURES}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"Atlas is missing columns: {missing}")
    analysis = data[["gene", TARGET, *FEATURES]].dropna().copy()
    if len(analysis) != len(data):
        raise ValueError(
            f"Technical analysis lost {len(data) - len(analysis)} genes to missing data"
        )

    predictions = cross_validate(analysis, args.folds, args.repeats, args.seed)
    repeats = repeat_metrics(predictions)
    averaged, bootstrap = bootstrap_average_predictions(
        predictions, args.bootstrap, args.seed
    )
    intervals = bootstrap_intervals(bootstrap)
    coefficients, full_alpha = full_model_coefficients(analysis)

    tables = ROOT / "outputs" / "tables"
    reports = ROOT / "outputs" / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(
        tables / "gene_predictability_technical_cv_predictions.csv.gz", index=False
    )
    repeats.to_csv(
        tables / "gene_predictability_technical_repeat_metrics.csv", index=False
    )
    averaged.to_csv(
        tables / "gene_predictability_technical_average_predictions.csv", index=False
    )
    bootstrap.to_csv(
        tables / "gene_predictability_technical_bootstrap.csv.gz", index=False
    )
    intervals.to_csv(
        tables / "gene_predictability_technical_bootstrap_intervals.csv", index=False
    )
    coefficients.to_csv(
        tables / "gene_predictability_technical_coefficients.csv", index=False
    )

    technical = intervals.loc[intervals["model"] == "technical"].set_index("metric")
    report = {
        "analysis": "held-out-gene technical predictability baseline",
        "target": TARGET,
        "genes": len(analysis),
        "features": FEATURES,
        "folds": args.folds,
        "repeats": args.repeats,
        "bootstrap_replicates": args.bootstrap,
        "full_model_alpha": full_alpha,
        "technical_model_intervals": technical.reset_index().to_dict(orient="records"),
        "repeat_metric_medians": (
            repeats.groupby("model")[["r2", "rmse", "mae", "spearman"]]
            .median()
            .reset_index()
            .to_dict(orient="records")
        ),
    }
    (reports / "gene_predictability_technical_model.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
