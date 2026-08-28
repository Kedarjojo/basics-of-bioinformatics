#!/usr/bin/env python3
"""Grouped ablation and permutation sensitivity analyses for UniProt biology."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import RepeatedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import model_gene_predictability_biology as core


ROOT = Path(__file__).resolve().parents[1]
LENGTH = ["log1p_protein_length"]
TOPOLOGY = ["log1p_transmembrane_count", "has_signal_peptide"]
LOCALIZATION = list(core.BINARY_BIOLOGICAL_FEATURES[1:])
# Remove broad terms that overlap strongly with more specific annotations.
REDUCED_LOCALIZATION = [
    feature for feature in LOCALIZATION
    if feature not in {"uniprot_membrane_any", "uniprot_extracellular"}
]

VARIANT_FEATURES = {
    "technical": core.TECHNICAL_FEATURES,
    "technical_plus_length": core.TECHNICAL_FEATURES + LENGTH,
    "technical_plus_topology": core.TECHNICAL_FEATURES + TOPOLOGY,
    "technical_plus_localization": core.TECHNICAL_FEATURES + LOCALIZATION,
    "combined_full": core.TECHNICAL_FEATURES + LENGTH + TOPOLOGY + LOCALIZATION,
    "combined_no_length": core.TECHNICAL_FEATURES + TOPOLOGY + LOCALIZATION,
    "combined_no_topology": core.TECHNICAL_FEATURES + LENGTH + LOCALIZATION,
    "combined_no_localization": core.TECHNICAL_FEATURES + LENGTH + TOPOLOGY,
    "combined_reduced_overlap": (
        core.TECHNICAL_FEATURES + LENGTH + TOPOLOGY + REDUCED_LOCALIZATION
    ),
}

COMPARISONS = [
    ("full_vs_technical", "combined_full", "technical"),
    ("reduced_vs_technical", "combined_reduced_overlap", "technical"),
    ("length_group_vs_technical", "technical_plus_length", "technical"),
    ("topology_group_vs_technical", "technical_plus_topology", "technical"),
    ("localization_group_vs_technical", "technical_plus_localization", "technical"),
    ("unique_length_full_vs_no_length", "combined_full", "combined_no_length"),
    ("unique_topology_full_vs_no_topology", "combined_full", "combined_no_topology"),
    (
        "unique_localization_full_vs_no_localization",
        "combined_full", "combined_no_localization",
    ),
    ("full_vs_reduced_overlap", "combined_full", "combined_reduced_overlap"),
]


def folds(n: int, n_splits: int, repeats: int, seed: int):
    splitter = RepeatedKFold(
        n_splits=n_splits, n_repeats=repeats, random_state=seed
    )
    return list(splitter.split(np.arange(n)))


def fit_variants(
    data: pd.DataFrame, split_list: list, variants: dict[str, list[str]]
) -> pd.DataFrame:
    y = data[core.TARGET].to_numpy(float)
    genes = data["gene"].astype(str).to_numpy()
    matrices = {
        name: data[features].to_numpy(float) for name, features in variants.items()
    }
    rows = []
    for split_index, (train, test) in enumerate(split_list):
        predictions = {}
        alphas = {}
        for name, matrix in matrices.items():
            model = make_pipeline(
                StandardScaler(), core.RidgeCV(alphas=core.ALPHAS)
            ).fit(matrix[train], y[train])
            predictions[name] = model.predict(matrix[test])
            alphas[name] = float(model[-1].alpha_)
        for local, position in enumerate(test):
            row = {"gene": genes[position], "split": split_index,
                   "observed": float(y[position])}
            for name in variants:
                row[f"pred_{name}"] = float(predictions[name][local])
                row[f"alpha_{name}"] = alphas[name]
            rows.append(row)
    return pd.DataFrame(rows)


def average_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    prediction_columns = [column for column in predictions if column.startswith("pred_")]
    aggregations = {"observed": ("observed", "first")}
    aggregations.update({column: (column, "mean") for column in prediction_columns})
    return predictions.groupby("gene", sort=True).agg(**aggregations).reset_index()


def bootstrap_comparisons(
    averaged: pd.DataFrame, replicates: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    metric_rows, delta_rows = [], []
    n = len(averaged)
    for replicate in range(replicates):
        draw = averaged.iloc[rng.integers(0, n, n)]
        observed = draw["observed"].to_numpy(float)
        values = {}
        for variant in VARIANT_FEATURES:
            values[variant] = core.metrics(
                observed, draw[f"pred_{variant}"].to_numpy(float)
            )
            metric_rows.append(
                {"replicate": replicate, "variant": variant, **values[variant]}
            )
        for comparison, enhanced, reference in COMPARISONS:
            for metric in ("r2", "rmse", "mae", "spearman"):
                delta = (
                    values[enhanced][metric] - values[reference][metric]
                    if metric in {"r2", "spearman"}
                    else values[reference][metric] - values[enhanced][metric]
                )
                delta_rows.append({
                    "replicate": replicate, "comparison": comparison,
                    "enhanced": enhanced, "reference": reference,
                    "metric": metric, "favorable_delta": float(delta),
                })
    return pd.DataFrame(metric_rows), pd.DataFrame(delta_rows)


def summarize_metrics(bootstrap: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for variant, group in bootstrap.groupby("variant", sort=True):
        for metric in ("r2", "rmse", "mae", "spearman"):
            values = group[metric].dropna().to_numpy(float)
            rows.append({
                "variant": variant, "metric": metric,
                "estimate": float(np.median(values)),
                "ci_low": float(np.quantile(values, 0.025)),
                "ci_high": float(np.quantile(values, 0.975)),
            })
    return pd.DataFrame(rows)


def summarize_deltas(deltas: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["comparison", "enhanced", "reference", "metric"]
    for key, group in deltas.groupby(keys, sort=True):
        values = group["favorable_delta"].to_numpy(float)
        rows.append(dict(zip(keys, key)) | {
            "estimate": float(np.median(values)),
            "ci_low": float(np.quantile(values, 0.025)),
            "ci_high": float(np.quantile(values, 0.975)),
            "probability_favorable": float(np.mean(values > 0)),
        })
    return pd.DataFrame(rows)


def averaged_metric(averaged: pd.DataFrame, variant: str) -> dict[str, float]:
    return core.metrics(
        averaged["observed"].to_numpy(float),
        averaged[f"pred_{variant}"].to_numpy(float),
    )


def permutation_test(
    data: pd.DataFrame,
    split_list: list,
    observed_average: pd.DataFrame,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    """Permute all biological rows jointly, preserving annotation correlation."""
    rng = np.random.default_rng(seed + 2)
    y = data[core.TARGET].to_numpy(float)
    technical = data[core.TECHNICAL_FEATURES].to_numpy(float)
    biological = data[core.BIOLOGICAL_FEATURES].to_numpy(float)
    technical_observed = averaged_metric(observed_average, "technical")
    full_observed = averaged_metric(observed_average, "combined_full")
    observed_delta = {
        metric: (
            full_observed[metric] - technical_observed[metric]
            if metric in {"r2", "spearman"}
            else technical_observed[metric] - full_observed[metric]
        )
        for metric in ("r2", "rmse", "mae", "spearman")
    }

    rows = []
    for permutation in range(permutations):
        shuffled = biological[rng.permutation(len(data))]
        matrix = np.column_stack([technical, shuffled])
        sums = np.zeros(len(data), dtype=float)
        counts = np.zeros(len(data), dtype=int)
        for train, test in split_list:
            model = make_pipeline(
                StandardScaler(), core.RidgeCV(alphas=core.ALPHAS)
            ).fit(matrix[train], y[train])
            sums[test] += model.predict(matrix[test])
            counts[test] += 1
        if np.any(counts == 0):
            raise RuntimeError("Permutation folds did not predict every gene")
        null_metrics = core.metrics(y, sums / counts)
        for metric in observed_delta:
            null_delta = (
                null_metrics[metric] - technical_observed[metric]
                if metric in {"r2", "spearman"}
                else technical_observed[metric] - null_metrics[metric]
            )
            rows.append({
                "permutation": permutation, "metric": metric,
                "null_favorable_delta": float(null_delta),
            })
    null = pd.DataFrame(rows)
    summary = {}
    for metric, observed in observed_delta.items():
        values = null.loc[null["metric"] == metric, "null_favorable_delta"].to_numpy()
        summary[metric] = {
            "observed_favorable_delta": float(observed),
            "null_median": float(np.median(values)),
            "null_95_high": float(np.quantile(values, 0.95)),
            "empirical_p_one_sided": float((1 + np.sum(values >= observed)) / (1 + len(values))),
        }
    return null, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, default=core.ATLAS)
    parser.add_argument("--annotations", type=Path, default=core.ANNOTATIONS)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260828)
    args = parser.parse_args()
    if min(args.folds, args.repeats, args.bootstrap, args.permutations) < 1:
        raise ValueError("All analysis counts must be positive")
    if args.folds < 2:
        raise ValueError("At least two folds are required")

    data = core.load_analysis(args.atlas, args.annotations)
    split_list = folds(len(data), args.folds, args.repeats, args.seed)
    predictions = fit_variants(data, split_list, VARIANT_FEATURES)
    averaged = average_predictions(predictions)
    bootstrap, deltas = bootstrap_comparisons(averaged, args.bootstrap, args.seed)
    model_summary = summarize_metrics(bootstrap)
    delta_summary = summarize_deltas(deltas)
    null, permutation_summary = permutation_test(
        data, split_list, averaged, args.permutations, args.seed
    )

    tables = ROOT / "outputs" / "tables"
    reports = ROOT / "outputs" / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    outputs = {
        "gene_predictability_biology_sensitivity_predictions.csv.gz": predictions,
        "gene_predictability_biology_sensitivity_average_predictions.csv": averaged,
        "gene_predictability_biology_sensitivity_bootstrap.csv.gz": bootstrap,
        "gene_predictability_biology_sensitivity_model_intervals.csv": model_summary,
        "gene_predictability_biology_sensitivity_paired_deltas.csv.gz": deltas,
        "gene_predictability_biology_sensitivity_paired_intervals.csv": delta_summary,
        "gene_predictability_biology_permutation_null.csv.gz": null,
    }
    for filename, frame in outputs.items():
        frame.to_csv(tables / filename, index=False)

    report = {
        "analysis": "biological predictability grouped ablation and sensitivity",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "target": core.TARGET, "genes": len(data),
        "folds": args.folds, "repeats": args.repeats,
        "bootstrap_replicates": args.bootstrap,
        "permutations": args.permutations, "seed": args.seed,
        "groups": {"length": LENGTH, "topology": TOPOLOGY,
                   "localization": LOCALIZATION,
                   "reduced_localization": REDUCED_LOCALIZATION},
        "variant_features": VARIANT_FEATURES,
        "comparisons": [
            {"name": name, "enhanced": enhanced, "reference": reference}
            for name, enhanced, reference in COMPARISONS
        ],
        "paired_intervals": delta_summary.to_dict(orient="records"),
        "permutation_test": permutation_summary,
        "interpretation": (
            "Positive favorable deltas favor the enhanced model. Biological rows are "
            "permuted jointly to preserve correlations among annotations while breaking "
            "their gene-specific relationship to the outcome and technical covariates."
        ),
    }
    (reports / "gene_predictability_biology_sensitivity.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
