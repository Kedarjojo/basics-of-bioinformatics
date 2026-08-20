"""Matched-estimator leave-one-cancer-out evaluation with auditable outputs."""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import REPORTS_DIR, TABLES_DIR
from .data_prep import load_model_table

RNA = ["tumor_log2tpm_mean", "rna_log2fc"]
GENE = ["gene_length", "protein_coding"]


def feature_sets(frame: pd.DataFrame) -> dict[str, list[str]]:
    localization = sorted(c for c in frame if c.startswith("loc_"))
    annotation = ["localization_available"]
    return {
        "rna": RNA,
        "rna_gene": RNA + GENE,
        "rna_localization": RNA + localization + annotation,
        "full": RNA + GENE + localization + annotation,
    }


def estimator(kind: str, features: list[str], seed: int) -> Pipeline:
    prep = ColumnTransformer(
        [("numeric", Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
        ]), features)],
        remainder="drop",
    )
    if kind == "logistic":
        model = LogisticRegression(
            class_weight="balanced", max_iter=2_000, random_state=seed
        )
    else:
        model = RandomForestClassifier(
            n_estimators=400,
            min_samples_leaf=5,
            class_weight="balanced",
            n_jobs=-1,
            random_state=seed,
        )
    return Pipeline([("prepare", prep), ("model", model)])


def metrics(y: pd.Series, probability: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "pr_auc": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def clustered_bootstrap(
    predictions: pd.DataFrame, comparisons: list[tuple[str, str]], n: int, seed: int
) -> pd.DataFrame:
    """Bootstrap paired metric differences by gene clusters."""
    if n <= 0:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    gene_indices = {
        gene: group.index.to_numpy()
        for gene, group in predictions.groupby("ensembl_base", sort=False)
    }
    genes = np.array(list(gene_indices), dtype=object)
    rows: list[dict] = []
    for replicate in range(n):
        sampled = rng.choice(genes, size=len(genes), replace=True)
        indices = np.concatenate([gene_indices[gene] for gene in sampled])
        boot = predictions.loc[indices]
        y = boot["protein_detected"]
        for reference, candidate in comparisons:
            for metric_name, scorer in (
                ("roc_auc", roc_auc_score),
                ("pr_auc", average_precision_score),
            ):
                delta = scorer(y, boot[f"prob_{candidate}"]) - scorer(
                    y, boot[f"prob_{reference}"]
                )
                rows.append(
                    {
                        "replicate": replicate,
                        "reference": reference,
                        "candidate": candidate,
                        "metric": metric_name,
                        "delta": float(delta),
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["logistic", "random_forest"], default="logistic")
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument(
        "--annotated-only",
        action="store_true",
        help="Restrict evaluation to genes with a localization annotation.",
    )
    args = parser.parse_args()

    frame, _ = load_model_table()
    if args.annotated_only:
        frame = frame[frame["localization_available"] == 1].copy()
    frame = frame.dropna(subset=["protein_detected", *RNA]).reset_index(drop=True)
    frame["protein_detected"] = frame["protein_detected"].astype(int)
    sets = feature_sets(frame)
    prediction_parts: list[pd.DataFrame] = []
    metric_rows: list[dict] = []

    for held_out in sorted(frame["cancer"].unique()):
        train = frame["cancer"] != held_out
        test = ~train
        fold_predictions = frame.loc[test, ["ensembl_base", "gene", "cancer", "protein_detected"]].copy()
        for name, features in sets.items():
            fitted = estimator(args.model, features, args.seed)
            fitted.fit(frame.loc[train, features], frame.loc[train, "protein_detected"])
            probability = fitted.predict_proba(frame.loc[test, features])[:, 1]
            fold_predictions[f"prob_{name}"] = probability
            metric_rows.append(
                {
                    "held_out_cancer": held_out,
                    "feature_set": name,
                    "model": args.model,
                    "annotated_only": args.annotated_only,
                    "n_test": int(test.sum()),
                    "prevalence": float(frame.loc[test, "protein_detected"].mean()),
                    **metrics(frame.loc[test, "protein_detected"], probability),
                }
            )
        prediction_parts.append(fold_predictions)
        print(f"Completed held-out cancer: {held_out}")

    predictions = pd.concat(prediction_parts, ignore_index=True)
    metric_table = pd.DataFrame(metric_rows)
    suffix = f"{args.model}_{'annotated' if args.annotated_only else 'all'}"
    predictions.to_csv(TABLES_DIR / f"loco_predictions_{suffix}.csv", index=False)
    metric_table.to_csv(TABLES_DIR / f"loco_metrics_{suffix}.csv", index=False)

    comparisons = [("rna", "rna_localization"), ("rna_gene", "full")]
    bootstrap = clustered_bootstrap(predictions, comparisons, args.bootstrap, args.seed)
    bootstrap.to_csv(TABLES_DIR / f"loco_bootstrap_{suffix}.csv", index=False)

    if len(bootstrap):
        bootstrap_summary = (
            bootstrap.groupby(["reference", "candidate", "metric"])["delta"]
            .agg(
                mean="mean",
                ci_low=lambda x: x.quantile(0.025),
                ci_high=lambda x: x.quantile(0.975),
                p_nonpositive=lambda x: (x <= 0).mean(),
            )
            .reset_index()
        )
    else:
        bootstrap_summary = pd.DataFrame()
    bootstrap_summary.to_csv(
        TABLES_DIR / f"loco_bootstrap_summary_{suffix}.csv", index=False
    )

    manifest = {
        "model": args.model,
        "seed": args.seed,
        "bootstrap_replicates": args.bootstrap,
        "annotated_only": args.annotated_only,
        "rows": len(frame),
        "ensembl_genes": int(frame["ensembl_base"].nunique()),
        "feature_sets": sets,
    }
    with (REPORTS_DIR / f"loco_manifest_{suffix}.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    print("\nMean metrics across held-out cancers:")
    print(metric_table.groupby("feature_set")[["roc_auc", "pr_auc", "brier"]].mean())
    if len(bootstrap_summary):
        print("\nGene-clustered paired bootstrap differences:")
        print(bootstrap_summary.to_string(index=False))


if __name__ == "__main__":
    main()

