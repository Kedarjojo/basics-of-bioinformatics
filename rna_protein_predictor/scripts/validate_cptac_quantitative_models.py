#!/usr/bin/env python3
"""Regularized, robust and bootstrap validation of the CPTAC benchmark."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import HuberRegressor, RidgeCV
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

try:
    from benchmark_cptac_quantitative import DEFAULT_CANCERS, load_cancer
except ModuleNotFoundError:  # Supports module-style execution in development.
    from src.benchmark_cptac_quantitative import DEFAULT_CANCERS, load_cancer


ROOT = Path(__file__).resolve().parents[1]
RIDGE_ALPHAS = np.array([0.01, 0.1, 1.0, 10.0, 100.0])


def _fit_gene(
    cancer: str,
    gene: str,
    x: np.ndarray,
    y: np.ndarray,
    sample_ids: np.ndarray,
    fold_ids: np.ndarray,
    min_pairs: int,
) -> list[dict]:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < min_pairs or np.ptp(y[valid]) == 0:
        return []
    rows: list[dict] = []
    for fold in sorted(np.unique(fold_ids)):
        train = valid & (fold_ids != fold)
        test = valid & (fold_ids == fold)
        if train.sum() < max(10, min_pairs // 2) or test.sum() == 0:
            continue
        x_train, y_train = x[train, None], y[train]
        x_test = x[test, None]
        baseline = float(np.mean(y_train))
        if np.ptp(x_train) == 0:
            pred_ols = np.full(test.sum(), baseline)
            pred_ridge = pred_ols.copy()
            pred_huber = pred_ols.copy()
            ridge_alpha = np.nan
        else:
            slope, intercept = np.polyfit(x_train[:, 0], y_train, 1)
            pred_ols = intercept + slope * x_test[:, 0]
            ridge = make_pipeline(
                StandardScaler(),
                RidgeCV(alphas=RIDGE_ALPHAS),
            ).fit(x_train, y_train)
            pred_ridge = ridge.predict(x_test)
            ridge_alpha = float(ridge[-1].alpha_)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=ConvergenceWarning)
                try:
                    huber = make_pipeline(
                        StandardScaler(),
                        HuberRegressor(epsilon=1.35, alpha=0.0001, max_iter=200),
                    ).fit(x_train, y_train)
                    pred_huber = huber.predict(x_test)
                except (ValueError, FloatingPointError):
                    pred_huber = pred_ridge.copy()
        test_positions = np.flatnonzero(test)
        train_min, train_max = float(np.min(y_train)), float(np.max(y_train))
        for j, position in enumerate(test_positions):
            rows.append(
                {
                    "cancer": cancer.upper(),
                    "gene": gene,
                    "sample_id": sample_ids[position],
                    "fold": int(fold),
                    "observed": float(y[position]),
                    "pred_gene_mean": baseline,
                    "pred_ols": float(pred_ols[j]),
                    "pred_ridge": float(pred_ridge[j]),
                    "pred_huber": float(pred_huber[j]),
                    "ridge_alpha": ridge_alpha,
                    "train_protein_min": train_min,
                    "train_protein_max": train_max,
                }
            )
    return rows


def fit_cancer(
    cancer: str,
    rna: pd.DataFrame,
    protein: pd.DataFrame,
    folds: int,
    seed: int,
    min_pairs: int,
    n_jobs: int,
) -> pd.DataFrame:
    fold_ids = np.zeros(len(rna), dtype=int)
    splitter = KFold(n_splits=folds, shuffle=True, random_state=seed)
    for fold, (_, test) in enumerate(splitter.split(np.arange(len(rna))), start=1):
        fold_ids[test] = fold
    sample_ids = rna.index.astype(str).to_numpy()
    results = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(_fit_gene)(
            cancer,
            gene,
            rna[gene].to_numpy(float),
            protein[gene].to_numpy(float),
            sample_ids,
            fold_ids,
            min_pairs,
        )
        for gene in rna.columns
    )
    return pd.DataFrame([row for gene_rows in results for row in gene_rows])


def add_errors(pred: pd.DataFrame) -> pd.DataFrame:
    for model in ("gene_mean", "ols", "ridge", "huber"):
        pred[f"se_{model}"] = (pred["observed"] - pred[f"pred_{model}"]) ** 2
        pred[f"ae_{model}"] = (pred["observed"] - pred[f"pred_{model}"]).abs()
        pred[f"outside_train_{model}"] = (
            (pred[f"pred_{model}"] < pred["train_protein_min"])
            | (pred[f"pred_{model}"] > pred["train_protein_max"])
        )
    return pred


def point_estimates(pred: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cancer, group in pred.groupby("cancer", sort=True):
        baseline_mse = group["se_gene_mean"].mean()
        for model in ("ols", "ridge", "huber"):
            rows.append(
                {
                    "cancer": cancer,
                    "model": model,
                    "rows": len(group),
                    "genes": group["gene"].nunique(),
                    "rmse": float(np.sqrt(group[f"se_{model}"].mean())),
                    "mae": float(group[f"ae_{model}"].mean()),
                    "relative_mse_improvement": float(
                        1 - group[f"se_{model}"].mean() / baseline_mse
                    ),
                    "outside_training_protein_range": float(
                        group[f"outside_train_{model}"].mean()
                    ),
                }
            )
    return pd.DataFrame(rows)


def patient_bootstrap(
    pred: pd.DataFrame, replicates: int, seed: int
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    error_cols = [f"se_{m}" for m in ("gene_mean", "ols", "ridge", "huber")]
    per_patient = pred.groupby(["cancer", "sample_id"])[error_cols].sum()
    rows = []
    for cancer, group in per_patient.groupby(level="cancer", sort=True):
        values = group.to_numpy(float)
        n = len(values)
        for replicate in range(replicates):
            sampled = values[rng.integers(0, n, n)].sum(axis=0)
            baseline = sampled[0]
            for j, model in enumerate(("ols", "ridge", "huber"), start=1):
                rows.append(
                    {
                        "cancer": cancer,
                        "replicate": replicate,
                        "model": model,
                        "relative_mse_improvement": 1 - sampled[j] / baseline,
                    }
                )
    return pd.DataFrame(rows)


def gene_bootstrap(pred: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed + 1)
    error_cols = [f"se_{m}" for m in ("gene_mean", "ols", "ridge", "huber")]
    per_gene = pred.groupby(["cancer", "gene"])[error_cols].sum().reset_index()
    rows = []
    for cancer, group in per_gene.groupby("cancer", sort=True):
        baseline = group["se_gene_mean"].to_numpy(float)
        improvements = {
            model: 1 - group[f"se_{model}"].to_numpy(float) / baseline
            for model in ("ols", "ridge", "huber")
        }
        n = len(group)
        for replicate in range(replicates):
            sampled = rng.integers(0, n, n)
            for model, values in improvements.items():
                draw = values[sampled]
                rows.append(
                    {
                        "cancer": cancer,
                        "replicate": replicate,
                        "model": model,
                        "median_gene_improvement": float(np.nanmedian(draw)),
                        "fraction_genes_improved": float(np.nanmean(draw > 0)),
                    }
                )
    return pd.DataFrame(rows)


def interval_summary(
    samples: pd.DataFrame, value_columns: list[str], kind: str
) -> pd.DataFrame:
    rows = []
    for (cancer, model), group in samples.groupby(["cancer", "model"], sort=True):
        for column in value_columns:
            values = group[column].dropna().to_numpy(float)
            rows.append(
                {
                    "bootstrap": kind,
                    "cancer": cancer,
                    "model": model,
                    "estimand": column,
                    "estimate": float(np.median(values)),
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancers", nargs="+", default=list(DEFAULT_CANCERS))
    parser.add_argument("--min-pairs", type=int, default=30)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260827)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--n-jobs", type=int, default=-1)
    args = parser.parse_args()

    predictions = []
    for cancer in args.cancers:
        rna, protein = load_cancer(cancer)
        predictions.append(
            fit_cancer(
                cancer, rna, protein, args.folds, args.seed,
                args.min_pairs, args.n_jobs,
            )
        )
    pred = add_errors(pd.concat(predictions, ignore_index=True))
    points = point_estimates(pred)
    patient = patient_bootstrap(pred, args.bootstrap, args.seed)
    gene = gene_bootstrap(pred, args.bootstrap, args.seed)
    intervals = pd.concat(
        [
            interval_summary(
                patient, ["relative_mse_improvement"], "patient"
            ),
            interval_summary(
                gene,
                ["median_gene_improvement", "fraction_genes_improved"],
                "gene",
            ),
        ],
        ignore_index=True,
    )

    tables = ROOT / "outputs" / "tables"
    reports = ROOT / "outputs" / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    pred.to_csv(tables / "cptac_quantitative_model_predictions.csv.gz", index=False)
    points.to_csv(tables / "cptac_quantitative_model_summary.csv", index=False)
    patient.to_csv(tables / "cptac_quantitative_patient_bootstrap.csv.gz", index=False)
    gene.to_csv(tables / "cptac_quantitative_gene_bootstrap.csv.gz", index=False)
    intervals.to_csv(tables / "cptac_quantitative_bootstrap_intervals.csv", index=False)

    overview = {
        "analysis": "regularized and robust CPTAC quantitative validation",
        "primary_model": "ols",
        "sensitivity_models": ["ridge", "huber"],
        "ridge_alphas": RIDGE_ALPHAS.tolist(),
        "bootstrap_replicates": args.bootstrap,
        "point_estimates": points.to_dict(orient="records"),
        "bootstrap_intervals": intervals.to_dict(orient="records"),
    }
    (reports / "cptac_quantitative_model_validation.json").write_text(
        json.dumps(overview, indent=2), encoding="utf-8"
    )
    print(json.dumps(overview, indent=2))


if __name__ == "__main__":
    main()
