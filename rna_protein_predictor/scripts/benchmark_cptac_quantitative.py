"""Prespecified quantitative RNA-to-protein benchmark using matched CPTAC data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import KFold


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CANCERS = ("brca", "coad", "gbm", "pdac")


def _safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    result = spearmanr(x, y) if method == "spearman" else pearsonr(x, y)
    return float(result.statistic)


def _linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 2 or np.ptp(x) == 0:
        return float(np.mean(y)), 0.0
    slope, intercept = np.polyfit(x, y, 1)
    return float(intercept), float(slope)


def load_cancer(cancer: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = ROOT / "data" / "cptac" / cancer
    rna_path = base / "matched_transcriptomics.csv.gz"
    protein_path = base / "matched_proteomics.csv.gz"
    sample_path = base / "matched_samples.csv"
    missing = [p for p in (rna_path, protein_path, sample_path) if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing CPTAC input(s): " + ", ".join(map(str, missing)))

    rna = pd.read_csv(rna_path, index_col=0)
    protein = pd.read_csv(protein_path, index_col=0)
    samples = pd.read_csv(sample_path)
    samples["sample_id"] = samples["sample_id"].astype(str)
    tumour_ids = samples.loc[
        samples["sample_kind_heuristic"] == "tumor_or_unspecified", "sample_id"
    ]
    common_samples = rna.index.astype(str).intersection(protein.index.astype(str))
    common_samples = common_samples.intersection(pd.Index(tumour_ids))
    common_genes = rna.columns.intersection(protein.columns)
    rna = rna.loc[common_samples, common_genes].apply(
    pd.to_numeric,
    errors="coerce",)
    protein = protein.loc[common_samples, common_genes].apply(
        pd.to_numeric,
        errors="coerce",)
    if (rna < 0).any().any():
        raise ValueError(
            "Negative WashU FPKM values found; log2(FPKM + 1) is not valid.")
    rna = np.log2(rna + 1.0)
    return rna, protein


def within_cancer_cv(
    cancer: str,
    rna: pd.DataFrame,
    protein: pd.DataFrame,
    min_pairs: int,
    folds: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    split = KFold(n_splits=folds, shuffle=True, random_state=seed)
    fold_ids = np.zeros(len(rna), dtype=int)
    for fold, (_, test) in enumerate(split.split(np.arange(len(rna))), start=1):
        fold_ids[test] = fold
    records: list[dict] = []
    for gene in rna.columns:
        x = rna[gene].to_numpy(float)
        y = protein[gene].to_numpy(float)
        valid = np.isfinite(x) & np.isfinite(y)
        if valid.sum() < min_pairs or np.ptp(y[valid]) == 0:
            continue
        for fold in range(1, folds + 1):
            train = np.flatnonzero(valid & (fold_ids != fold))
            test = np.flatnonzero(valid & (fold_ids == fold))
            if len(train) < 2 or len(test) == 0:
                continue
            intercept, slope = _linear_fit(x[train], y[train])
            baseline = float(np.mean(y[train]))
            for idx in test:
                records.append(
                    {
                        "cancer": cancer.upper(),
                        "gene": gene,
                        "sample_id": rna.index[idx],
                        "fold": fold,
                        "observed": y[idx],
                        "pred_gene_mean": baseline,
                        "pred_rna_linear": intercept + slope * x[idx],
                    }
                )
    pred = pd.DataFrame.from_records(records)
    metrics: list[dict] = []
    for (cancer_code, gene), group in pred.groupby(["cancer", "gene"], sort=True):
        observed = group["observed"].to_numpy()
        base = group["pred_gene_mean"].to_numpy()
        fitted = group["pred_rna_linear"].to_numpy()
        mse_base = float(np.mean((observed - base) ** 2))
        mse_rna = float(np.mean((observed - fitted) ** 2))
        metrics.append(
            {
                "cancer": cancer_code,
                "gene": gene,
                "n_pairs": len(group),
                "mse_gene_mean": mse_base,
                "mse_rna_linear": mse_rna,
                "relative_mse_improvement": 1 - mse_rna / mse_base if mse_base > 0 else np.nan,
                "spearman_observed_rna_prediction": _safe_corr(observed, fitted, "spearman"),
                "pearson_observed_rna_prediction": _safe_corr(observed, fitted, "pearson"),
            }
        )
    return pred, pd.DataFrame(metrics)


def summarize_cv(pred: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cancer, group in pred.groupby("cancer", sort=True):
        observed = group["observed"].to_numpy()
        base = group["pred_gene_mean"].to_numpy()
        fitted = group["pred_rna_linear"].to_numpy()
        mse_base = float(np.mean((observed - base) ** 2))
        mse_rna = float(np.mean((observed - fitted) ** 2))
        gm = genes.loc[genes["cancer"] == cancer]
        rows.append(
            {
                "cancer": cancer,
                "prediction_rows": len(group),
                "genes": int(gm["gene"].nunique()),
                "rmse_gene_mean": np.sqrt(mse_base),
                "rmse_rna_linear": np.sqrt(mse_rna),
                "relative_mse_improvement": 1 - mse_rna / mse_base,
                "median_gene_spearman": gm["spearman_observed_rna_prediction"].median(),
                "fraction_genes_positive_spearman": (gm["spearman_observed_rna_prediction"] > 0).mean(),
                "fraction_genes_positive_mse_improvement": (gm["relative_mse_improvement"] > 0).mean(),
            }
        )
    return pd.DataFrame(rows)


def observed_coupling(
    matrices: dict[str, tuple[pd.DataFrame, pd.DataFrame]], min_pairs: int
) -> pd.DataFrame:
    rows = []
    for cancer, (rna, protein) in matrices.items():
        for gene in rna.columns:
            x, y = rna[gene].to_numpy(float), protein[gene].to_numpy(float)
            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() < min_pairs:
                continue
            rows.append(
                {
                    "cancer": cancer.upper(),
                    "gene": gene,
                    "n_pairs": int(valid.sum()),
                    "spearman_rna_protein": _safe_corr(x[valid], y[valid], "spearman"),
                    "pearson_rna_protein": _safe_corr(x[valid], y[valid], "pearson"),
                }
            )
    return pd.DataFrame(rows)


def coupling_stability(coupling: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide = coupling.pivot(index="gene", columns="cancer", values="spearman_rna_protein")
    pairwise = wide.corr(method="spearman", min_periods=100).rename_axis(
        index="cancer_a", columns="cancer_b"
    )
    pairs = (
        pairwise.where(np.triu(np.ones(pairwise.shape), k=1).astype(bool))
        .stack()
        .rename("gene_profile_spearman")
        .reset_index()
    )
    gene = pd.DataFrame(
        {
            "gene": wide.index,
            "cancers_observed": wide.notna().sum(axis=1),
            "median_spearman": wide.median(axis=1),
            "min_spearman": wide.min(axis=1),
            "max_spearman": wide.max(axis=1),
            "positive_cancers": (wide > 0).sum(axis=1),
            "negative_cancers": (wide < 0).sum(axis=1),
        }
    ).reset_index(drop=True)
    return pairs, gene


def loco_transfer(
    matrices: dict[str, tuple[pd.DataFrame, pd.DataFrame]], min_pairs: int
) -> pd.DataFrame:
    rows = []
    for held, (test_rna, test_protein) in matrices.items():
        training = {k: v for k, v in matrices.items() if k != held}
        for gene in test_rna.columns:
            slopes = []
            for train_rna, train_protein in training.values():
                if gene not in train_rna or gene not in train_protein:
                    continue
                x, y = train_rna[gene].to_numpy(float), train_protein[gene].to_numpy(float)
                valid = np.isfinite(x) & np.isfinite(y)
                if valid.sum() < min_pairs or np.std(x[valid]) == 0 or np.std(y[valid]) == 0:
                    continue
                xz = (x[valid] - np.mean(x[valid])) / np.std(x[valid], ddof=0)
                yz = (y[valid] - np.mean(y[valid])) / np.std(y[valid], ddof=0)
                slopes.append(float(np.dot(xz, yz) / np.dot(xz, xz)))
            x, y = test_rna[gene].to_numpy(float), test_protein[gene].to_numpy(float)
            valid = np.isfinite(x) & np.isfinite(y)
            if len(slopes) < 2 or valid.sum() < min_pairs or np.std(x[valid]) == 0:
                continue
            predicted = np.mean(slopes) * (
                (x[valid] - np.mean(x[valid])) / np.std(x[valid], ddof=0)
            )
            rows.append(
                {
                    "held_out_cancer": held.upper(),
                    "gene": gene,
                    "n_test_pairs": int(valid.sum()),
                    "training_cancers": len(slopes),
                    "mean_training_slope": float(np.mean(slopes)),
                    "test_spearman": _safe_corr(predicted, y[valid], "spearman"),
                    "test_pearson": _safe_corr(predicted, y[valid], "pearson"),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancers", nargs="+", default=list(DEFAULT_CANCERS))
    parser.add_argument("--min-pairs", type=int, default=30)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260827)
    args = parser.parse_args()

    matrices = {c: load_cancer(c) for c in args.cancers}
    predictions, gene_metrics = [], []
    for cancer, (rna, protein) in matrices.items():
        pred, metrics = within_cancer_cv(
            cancer, rna, protein, args.min_pairs, args.folds, args.seed
        )
        predictions.append(pred)
        gene_metrics.append(metrics)
    predictions = pd.concat(predictions, ignore_index=True)
    gene_metrics = pd.concat(gene_metrics, ignore_index=True)
    summary = summarize_cv(predictions, gene_metrics)
    coupling = observed_coupling(matrices, args.min_pairs)
    pairwise, stability = coupling_stability(coupling)
    loco = loco_transfer(matrices, args.min_pairs)

    table_dir = ROOT / "outputs" / "tables"
    report_dir = ROOT / "outputs" / "reports"
    table_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(table_dir / "cptac_quantitative_cv_predictions.csv.gz", index=False)
    gene_metrics.to_csv(table_dir / "cptac_quantitative_gene_metrics.csv", index=False)
    summary.to_csv(table_dir / "cptac_quantitative_cv_summary.csv", index=False)
    coupling.to_csv(table_dir / "cptac_observed_gene_coupling.csv", index=False)
    pairwise.to_csv(table_dir / "cptac_coupling_pairwise_stability.csv", index=False)
    stability.to_csv(table_dir / "cptac_coupling_gene_stability.csv", index=False)
    loco.to_csv(table_dir / "cptac_quantitative_loco_transfer.csv", index=False)

    loco_summary = (
        loco.groupby("held_out_cancer")
        .agg(
            genes=("gene", "nunique"),
            median_test_spearman=("test_spearman", "median"),
            fraction_positive=("test_spearman", lambda x: float((x > 0).mean())),
        )
        .reset_index()
    )
    overview = {
        "analysis": "prespecified quantitative CPTAC benchmark",
        "cancers": args.cancers,
        "min_pairs": args.min_pairs,
        "folds": args.folds,
        "seed": args.seed,
        "within_cancer_cv": summary.to_dict(orient="records"),
        "loco_transfer": loco_summary.to_dict(orient="records"),
        "pairwise_gene_profile_stability": pairwise.to_dict(orient="records"),
    }
    (report_dir / "cptac_quantitative_benchmark.json").write_text(
        json.dumps(overview, indent=2), encoding="utf-8"
    )
    print(json.dumps(overview, indent=2))


if __name__ == "__main__":
    main()
