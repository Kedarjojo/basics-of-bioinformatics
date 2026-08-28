#!/usr/bin/env python3
"""Test biological explanations of CPTAC gene predictability.

Compares intercept-only, technical-only, biological-only, and combined Ridge
models using identical repeated held-out-gene folds. Biological annotations are
joined from the frozen synonym-aware UniProt audit. All model comparisons use
the same complete-case gene set and paired bootstrap draws.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
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
ANNOTATIONS = ROOT / "outputs" / "tables" / "cptac_atlas_uniprot_annotations.csv"
TARGET = "primary_median_ridge_improvement"
TECHNICAL_FEATURES = [
    "median_rna_mean", "median_rna_sd", "median_rna_iqr",
    "median_protein_mean", "median_protein_sd", "median_protein_iqr",
    "median_protein_coverage", "minimum_protein_coverage",
    "median_paired_patients", "cancers_observed",
]
BINARY_BIOLOGICAL_FEATURES = [
    "has_signal_peptide",
    "uniprot_secreted", "uniprot_extracellular",
    "uniprot_plasma_membrane", "uniprot_membrane_any",
    "uniprot_cytosol", "uniprot_nucleus", "uniprot_mitochondrion",
    "uniprot_endoplasmic_reticulum", "uniprot_golgi",
    "uniprot_lysosome", "uniprot_peroxisome", "uniprot_cell_junction",
]
DERIVED_BIOLOGICAL_FEATURES = ["log1p_protein_length", "log1p_transmembrane_count"]
BIOLOGICAL_FEATURES = DERIVED_BIOLOGICAL_FEATURES + BINARY_BIOLOGICAL_FEATURES
MODEL_FEATURES = {
    "technical": TECHNICAL_FEATURES,
    "biological": BIOLOGICAL_FEATURES,
    "combined": TECHNICAL_FEATURES + BIOLOGICAL_FEATURES,
}
MODELS = ["intercept_only", "technical", "biological", "combined"]
ALPHAS = np.logspace(-4, 4, 17)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(spearmanr(x, y).statistic)


def metrics(observed: np.ndarray, fitted: np.ndarray) -> dict[str, float]:
    return {
        "r2": float(r2_score(observed, fitted)),
        "rmse": float(np.sqrt(mean_squared_error(observed, fitted))),
        "mae": float(mean_absolute_error(observed, fitted)),
        "spearman": safe_spearman(observed, fitted),
    }


def coerce_boolean(series: pd.Series, column: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    normalized = series.astype(str).str.strip().str.lower()
    mapping = {"true": True, "false": False, "1": True, "0": False}
    converted = normalized.map(mapping)
    if converted.isna().any():
        bad = sorted(normalized.loc[converted.isna()].unique())[:10]
        raise ValueError(f"Cannot parse Boolean column {column}: {bad}")
    return converted.astype(bool)


def load_analysis(atlas_path: Path, annotation_path: Path) -> pd.DataFrame:
    atlas = pd.read_csv(atlas_path)
    annotations = pd.read_csv(annotation_path)
    required_atlas = {"gene", TARGET, *TECHNICAL_FEATURES}
    required_annotations = {
        "gene", "mapping_method", "uniprot_annotation_available",
        "protein_length", "transmembrane_domain_count",
        *BINARY_BIOLOGICAL_FEATURES,
    }
    missing_atlas = sorted(required_atlas.difference(atlas.columns))
    missing_annotations = sorted(required_annotations.difference(annotations.columns))
    if missing_atlas:
        raise ValueError(f"Atlas is missing columns: {missing_atlas}")
    if missing_annotations:
        raise ValueError(f"Annotation table is missing columns: {missing_annotations}")
    if atlas["gene"].duplicated().any() or annotations["gene"].duplicated().any():
        raise ValueError("Atlas and annotation tables must contain unique gene rows")

    selected_annotations = annotations[[
        "gene", "mapping_method", "uniprot_annotation_available",
        "protein_length", "transmembrane_domain_count",
        *BINARY_BIOLOGICAL_FEATURES,
    ]].copy()
    data = atlas[["gene", TARGET, *TECHNICAL_FEATURES]].merge(
        selected_annotations, on="gene", how="left", validate="one_to_one"
    )
    available_raw = data["uniprot_annotation_available"].fillna(False)
    available = coerce_boolean(available_raw, "uniprot_annotation_available")
    data = data.loc[available].copy()
    data["log1p_protein_length"] = np.log1p(data["protein_length"].astype(float))
    data["log1p_transmembrane_count"] = np.log1p(
        data["transmembrane_domain_count"].astype(float)
    )
    for feature in BINARY_BIOLOGICAL_FEATURES:
        data[feature] = coerce_boolean(data[feature], feature).astype(int)
    required_complete = [TARGET, *TECHNICAL_FEATURES, *BIOLOGICAL_FEATURES]
    incomplete = data[required_complete].isna().any(axis=1)
    if incomplete.any():
        genes = data.loc[incomplete, "gene"].head(20).tolist()
        raise ValueError(
            f"Annotated common set has {int(incomplete.sum())} incomplete rows; "
            f"examples: {genes}"
        )
    return data.reset_index(drop=True)


def cross_validate(
    data: pd.DataFrame, folds: int, repeats: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    y = data[TARGET].to_numpy(float)
    genes = data["gene"].astype(str).to_numpy()
    matrices = {
        model: data[features].to_numpy(float)
        for model, features in MODEL_FEATURES.items()
    }
    splitter = RepeatedKFold(
        n_splits=folds, n_repeats=repeats, random_state=seed
    )
    prediction_rows: list[dict] = []
    coefficient_rows: list[dict] = []
    for split_index, (train, test) in enumerate(splitter.split(y)):
        repeat, fold = divmod(split_index, folds)
        predictions = {"intercept_only": np.full(len(test), y[train].mean())}
        for model_name, features in MODEL_FEATURES.items():
            model = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))
            model.fit(matrices[model_name][train], y[train])
            predictions[model_name] = model.predict(matrices[model_name][test])
            alpha = float(model[-1].alpha_)
            for feature, coefficient in zip(features, model[-1].coef_):
                coefficient_rows.append({
                    "repeat": repeat, "fold": fold, "model": model_name,
                    "feature": feature,
                    "standardized_coefficient": float(coefficient),
                    "ridge_alpha": alpha,
                })
        for local_index, position in enumerate(test):
            row = {
                "gene": genes[position], "repeat": repeat, "fold": fold,
                "observed": float(y[position]),
            }
            for model_name in MODELS:
                row[f"pred_{model_name}"] = float(predictions[model_name][local_index])
            prediction_rows.append(row)
    return pd.DataFrame(prediction_rows), pd.DataFrame(coefficient_rows)


def repeat_metrics(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for repeat, group in predictions.groupby("repeat", sort=True):
        observed = group["observed"].to_numpy(float)
        for model in MODELS:
            rows.append({"repeat": int(repeat), "model": model,
                         **metrics(observed, group[f"pred_{model}"].to_numpy(float))})
    return pd.DataFrame(rows)


def average_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    aggregations = {"observed": ("observed", "first")}
    for model in MODELS:
        aggregations[f"pred_{model}"] = (f"pred_{model}", "mean")
    return predictions.groupby("gene", sort=True).agg(**aggregations).reset_index()


def bootstrap_models_and_deltas(
    averaged: pd.DataFrame, replicates: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    model_rows, delta_rows = [], []
    comparisons = [
        ("biological_vs_intercept", "biological", "intercept_only"),
        ("combined_vs_technical", "combined", "technical"),
        ("combined_vs_biological", "combined", "biological"),
    ]
    n = len(averaged)
    for replicate in range(replicates):
        draw = averaged.iloc[rng.integers(0, n, n)]
        observed = draw["observed"].to_numpy(float)
        values = {}
        for model in MODELS:
            values[model] = metrics(observed, draw[f"pred_{model}"].to_numpy(float))
            model_rows.append({"replicate": replicate, "model": model, **values[model]})
        for comparison, enhanced, reference in comparisons:
            for metric in ("r2", "rmse", "mae", "spearman"):
                # Positive always favors the enhanced model.
                delta = (
                    values[enhanced][metric] - values[reference][metric]
                    if metric in {"r2", "spearman"}
                    else values[reference][metric] - values[enhanced][metric]
                )
                delta_rows.append({
                    "replicate": replicate, "comparison": comparison,
                    "metric": metric, "favorable_delta": float(delta),
                })
    return pd.DataFrame(model_rows), pd.DataFrame(delta_rows)


def model_intervals(bootstrap: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, group in bootstrap.groupby("model", sort=True):
        for metric in ("r2", "rmse", "mae", "spearman"):
            values = group[metric].dropna().to_numpy(float)
            rows.append({
                "model": model, "metric": metric,
                "estimate": float(np.median(values)),
                "ci_low": float(np.quantile(values, 0.025)),
                "ci_high": float(np.quantile(values, 0.975)),
            })
    return pd.DataFrame(rows)


def delta_intervals(deltas: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (comparison, metric), group in deltas.groupby(["comparison", "metric"]):
        values = group["favorable_delta"].to_numpy(float)
        rows.append({
            "comparison": comparison, "metric": metric,
            "estimate": float(np.median(values)),
            "ci_low": float(np.quantile(values, 0.025)),
            "ci_high": float(np.quantile(values, 0.975)),
            "probability_favorable": float(np.mean(values > 0)),
        })
    return pd.DataFrame(rows)


def coefficient_stability(coefficients: pd.DataFrame) -> pd.DataFrame:
    return (
        coefficients.groupby(["model", "feature"], sort=True)
        .agg(
            median_standardized_coefficient=("standardized_coefficient", "median"),
            ci_low=("standardized_coefficient", lambda x: x.quantile(0.025)),
            ci_high=("standardized_coefficient", lambda x: x.quantile(0.975)),
            positive_fraction=("standardized_coefficient", lambda x: (x > 0).mean()),
            negative_fraction=("standardized_coefficient", lambda x: (x < 0).mean()),
            fits=("standardized_coefficient", "size"),
        )
        .reset_index()
    )


def full_coefficients(data: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float]]:
    rows, alphas = [], {}
    y = data[TARGET].to_numpy(float)
    for model_name, features in MODEL_FEATURES.items():
        model = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))
        model.fit(data[features].to_numpy(float), y)
        alphas[model_name] = float(model[-1].alpha_)
        for feature, coefficient in zip(features, model[-1].coef_):
            rows.append({"model": model_name, "feature": feature,
                         "standardized_coefficient": float(coefficient)})
    return pd.DataFrame(rows), alphas


def descriptive_groups(
    data: pd.DataFrame, replicates: int, seed: int
) -> pd.DataFrame:
    """Unadjusted group estimates for plotting; not causal model coefficients."""
    rng = np.random.default_rng(seed + 1)
    rows = []
    for feature in BINARY_BIOLOGICAL_FEATURES:
        absent = data.loc[data[feature] == 0, TARGET].to_numpy(float)
        present = data.loc[data[feature] == 1, TARGET].to_numpy(float)
        differences = np.empty(replicates)
        for i in range(replicates):
            a = rng.choice(absent, len(absent), replace=True)
            p = rng.choice(present, len(present), replace=True)
            differences[i] = np.median(p) - np.median(a)
        rows.append({
            "feature": feature, "absent_n": len(absent), "present_n": len(present),
            "absent_median": float(np.median(absent)),
            "present_median": float(np.median(present)),
            "present_minus_absent_median": float(np.median(present) - np.median(absent)),
            "bootstrap_ci_low": float(np.quantile(differences, 0.025)),
            "bootstrap_ci_high": float(np.quantile(differences, 0.975)),
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, default=ATLAS)
    parser.add_argument("--annotations", type=Path, default=ANNOTATIONS)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260828)
    args = parser.parse_args()
    if args.folds < 2 or args.repeats < 1 or args.bootstrap < 1:
        raise ValueError("folds >= 2, repeats >= 1, and bootstrap >= 1 are required")

    data = load_analysis(args.atlas, args.annotations)
    predictions, fold_coefficients = cross_validate(
        data, args.folds, args.repeats, args.seed
    )
    repeats = repeat_metrics(predictions)
    averaged = average_predictions(predictions)
    bootstrap, deltas = bootstrap_models_and_deltas(
        averaged, args.bootstrap, args.seed
    )
    intervals = model_intervals(bootstrap)
    delta_summary = delta_intervals(deltas)
    stability = coefficient_stability(fold_coefficients)
    coefficients, full_alphas = full_coefficients(data)
    groups = descriptive_groups(data, args.bootstrap, args.seed)

    tables = ROOT / "outputs" / "tables"
    reports = ROOT / "outputs" / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    outputs = {
        "gene_predictability_biology_cv_predictions.csv.gz": predictions,
        "gene_predictability_biology_repeat_metrics.csv": repeats,
        "gene_predictability_biology_average_predictions.csv": averaged,
        "gene_predictability_biology_bootstrap.csv.gz": bootstrap,
        "gene_predictability_biology_model_intervals.csv": intervals,
        "gene_predictability_biology_paired_bootstrap_deltas.csv.gz": deltas,
        "gene_predictability_biology_paired_intervals.csv": delta_summary,
        "gene_predictability_biology_fold_coefficients.csv.gz": fold_coefficients,
        "gene_predictability_biology_coefficient_stability.csv": stability,
        "gene_predictability_biology_full_coefficients.csv": coefficients,
        "gene_predictability_biology_group_summary.csv": groups,
    }
    for filename, frame in outputs.items():
        frame.to_csv(tables / filename, index=False)

    report = {
        "analysis": "held-out-gene technical and biological predictability comparison",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "target": TARGET,
        "genes": len(data),
        "mapping_method_counts": {
            str(k): int(v) for k, v in data["mapping_method"].value_counts().items()
        },
        "technical_features": TECHNICAL_FEATURES,
        "biological_features": BIOLOGICAL_FEATURES,
        "models": MODELS,
        "folds": args.folds, "repeats": args.repeats,
        "bootstrap_replicates": args.bootstrap, "seed": args.seed,
        "ridge_alphas": ALPHAS.tolist(), "full_model_alphas": full_alphas,
        "model_intervals": intervals.to_dict(orient="records"),
        "paired_model_intervals": delta_summary.to_dict(orient="records"),
        "interpretation_rule": (
            "Positive favorable_delta favors the first model named in comparison. "
            "Biological incremental value is assessed primarily by combined_vs_technical."
        ),
        "descriptive_group_warning": (
            "Localization group summaries are unadjusted descriptions; adjusted inference "
            "comes from held-out combined-versus-technical prediction."
        ),
    }
    (reports / "gene_predictability_biology_model.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
