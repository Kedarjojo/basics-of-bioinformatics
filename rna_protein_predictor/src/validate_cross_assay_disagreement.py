"""Threshold sensitivity and gene-clustered validation of cross-assay disagreement."""

from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.cov_struct import Exchangeable
from statsmodels.stats.multitest import multipletests

from .config import DATA_PROCESSED, REPORTS_DIR, TABLES_DIR


def sensitivity_grid(atlas: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for hpa_threshold in (0.0, 0.10, 0.25, 0.49):
        for coverage_threshold in (0.50, 0.80, 0.90, 1.00):
            working = atlas.copy()
            working["event"] = (
                working["detected_fraction"].le(hpa_threshold)
                & working["cptac_protein_coverage"].ge(coverage_threshold)
            )
            recurrence = working.groupby("gene")["event"].sum()
            for cancer, subset in working.groupby("cancer_code"):
                hpa_low = subset["detected_fraction"].le(hpa_threshold)
                rows.append({
                    "hpa_detection_fraction_max": hpa_threshold,
                    "cptac_coverage_min": coverage_threshold,
                    "cancer_code": cancer,
                    "overlap_genes": len(subset),
                    "hpa_low_genes": int(hpa_low.sum()),
                    "events": int(subset["event"].sum()),
                    "event_rate_among_overlap": float(subset["event"].mean()),
                    "event_rate_among_hpa_low": (
                        float(subset.loc[hpa_low, "event"].mean()) if hpa_low.any() else np.nan
                    ),
                    "genes_recurrent_2plus": int((recurrence >= 2).sum()),
                    "genes_recurrent_3plus": int((recurrence >= 3).sum()),
                    "genes_recurrent_4": int((recurrence >= 4).sum()),
                })
    return pd.DataFrame(rows)


def gee_localization(atlas: pd.DataFrame, loc_cols: list[str]) -> pd.DataFrame:
    annotated = atlas[atlas[loc_cols].fillna(0).sum(axis=1) > 0].copy()
    annotated["log_gene_length"] = np.log1p(annotated["gene_length"])
    annotated["log_hpa_patients"] = np.log1p(annotated["total_patients"])
    rows = []
    specifications = {
        "hpa_negative_among_cptac_quantified": {
            "subset": annotated[annotated["cptac_quantified"]].copy(),
            "outcome": "hpa_negative",
            "covariates": "C(cancer_code) + log_gene_length + log_hpa_patients + cptac_protein_coverage",
        },
        "cptac_quantified_among_hpa_negative": {
            "subset": annotated[annotated["hpa_negative"]].copy(),
            "outcome": "cptac_quantified",
            "covariates": "C(cancer_code) + log_gene_length + log_hpa_patients",
        },
    }
    for analysis, spec in specifications.items():
        frame = spec["subset"].dropna(
            subset=[spec["outcome"], "gene", "log_gene_length", "log_hpa_patients"]
        ).copy()
        frame[spec["outcome"]] = frame[spec["outcome"]].astype(int)
        for feature in loc_cols:
            if frame[feature].nunique() < 2:
                continue
            formula = f"{spec['outcome']} ~ {feature} + {spec['covariates']}"
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    model = smf.gee(
                        formula=formula,
                        groups="gene",
                        data=frame,
                        family=sm.families.Binomial(),
                        cov_struct=Exchangeable(),
                    ).fit(maxiter=200)
                coefficient = float(model.params[feature])
                standard_error = float(model.bse[feature])
                rows.append({
                    "analysis": analysis,
                    "localization": feature,
                    "n_gene_cancer_pairs": len(frame),
                    "n_gene_clusters": frame["gene"].nunique(),
                    "coefficient": coefficient,
                    "standard_error": standard_error,
                    "odds_ratio": float(np.exp(coefficient)),
                    "ci_low": float(np.exp(coefficient - 1.96 * standard_error)),
                    "ci_high": float(np.exp(coefficient + 1.96 * standard_error)),
                    "pvalue": float(model.pvalues[feature]),
                })
            except (ValueError, np.linalg.LinAlgError, KeyError) as error:
                rows.append({
                    "analysis": analysis,
                    "localization": feature,
                    "error": str(error),
                })
    results = pd.DataFrame(rows)
    valid = results["pvalue"].notna()
    results["fdr"] = np.nan
    if valid.any():
        results.loc[valid, "fdr"] = (
            results.loc[valid]
            .groupby("analysis")["pvalue"]
            .transform(lambda x: multipletests(x, method="fdr_bh")[1])
        )
    return results


def main() -> None:
    atlas = pd.read_csv(TABLES_DIR / "cross_cancer_hpa_cptac_atlas.csv")
    master = pd.read_csv(DATA_PROCESSED / "master_table_augmented.csv")
    gene_features = (
        master[["gene", "gene_length", "protein_coding"]]
        .dropna(subset=["gene"])
        .groupby("gene", as_index=False)
        .agg(gene_length=("gene_length", "median"), protein_coding=("protein_coding", "max"))
    )
    atlas = atlas.drop(columns=["gene_length", "protein_coding"], errors="ignore")
    atlas = atlas.merge(gene_features, on="gene", how="left", validate="many_to_one")
    loc_cols = [column for column in atlas if column.startswith("loc_")]

    sensitivity = sensitivity_grid(atlas)
    gee = gee_localization(atlas, loc_cols)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    sensitivity.to_csv(TABLES_DIR / "cross_assay_threshold_sensitivity.csv", index=False)
    gee.to_csv(TABLES_DIR / "cross_assay_localization_gee.csv", index=False)

    primary = sensitivity[
        (sensitivity["hpa_detection_fraction_max"] == 0.49)
        & (sensitivity["cptac_coverage_min"] == 0.80)
    ]
    significant = gee[gee["fdr"] < 0.05]
    overview = {
        "primary_threshold": {
            "hpa_detection_fraction_max": 0.49,
            "cptac_coverage_min": 0.80,
            "by_cancer": primary.set_index("cancer_code")[
                ["events", "event_rate_among_overlap", "event_rate_among_hpa_low"]
            ].to_dict(orient="index"),
        },
        "sensitivity_combinations": int(
            sensitivity[["hpa_detection_fraction_max", "cptac_coverage_min"]]
            .drop_duplicates().shape[0]
        ),
        "gee_annotated_genes_only": True,
        "gee_gene_clustered": True,
        "gee_significant_localization_effects": int(len(significant)),
    }
    (REPORTS_DIR / "cross_assay_validation_overview.json").write_text(
        json.dumps(overview, indent=2), encoding="utf-8"
    )
    print(json.dumps(overview, indent=2))


if __name__ == "__main__":
    main()

