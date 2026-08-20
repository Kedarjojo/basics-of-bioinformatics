"""Build the HPA-IHC versus CPTAC-MS cross-cancer disagreement atlas."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests

from .config import DATA_PROCESSED, REPORTS_DIR, TABLES_DIR

ROOT = Path(__file__).resolve().parents[1]
COHORTS = {"coad": "COAD", "brca": "BRCA", "gbm": "GBM", "pdac": "PAAD"}


def gene_metrics(rna: pd.DataFrame, protein: pd.DataFrame, samples: pd.Index) -> pd.DataFrame:
    genes = rna.columns.intersection(protein.columns)
    rows = []
    for gene in genes:
        x = pd.to_numeric(rna.loc[samples, gene], errors="coerce").to_numpy(float)
        y = pd.to_numeric(protein.loc[samples, gene], errors="coerce").to_numpy(float)
        complete = np.isfinite(x) & np.isfinite(y)
        n = int(complete.sum())
        rho = pvalue = np.nan
        if n >= 30 and np.std(x[complete]) > 0 and np.std(y[complete]) > 0:
            rho, pvalue = spearmanr(x[complete], y[complete])
        rows.append({
            "gene": str(gene).upper(),
            "n_tumors": len(samples),
            "n_protein_quantified": n,
            "cptac_protein_coverage": n / len(samples),
            "rna_protein_spearman": rho,
            "rna_protein_p": pvalue,
        })
    return pd.DataFrame(rows)


def fdr(values: pd.Series) -> np.ndarray:
    valid = values.notna()
    result = np.full(len(values), np.nan)
    if valid.any():
        result[valid] = multipletests(values[valid], method="fdr_bh")[1]
    return result


def localization_tests(atlas: pd.DataFrame, loc_cols: list[str]) -> pd.DataFrame:
    """Decompose disagreement from both assay perspectives."""
    rows = []
    groups = [("ALL", atlas), *list(atlas.groupby("cancer_code"))]
    for cancer, frame in groups:
        analyses = {
            # Conditional on MS quantification, test IHC-negative enrichment.
            "hpa_negative_among_cptac_quantified": frame[frame["cptac_quantified"]],
            # Conditional on HPA negativity, test MS-quantification enrichment.
            "cptac_quantified_among_hpa_negative": frame[frame["hpa_negative"]],
        }
        outcomes = {
            "hpa_negative_among_cptac_quantified": "hpa_negative",
            "cptac_quantified_among_hpa_negative": "cptac_quantified",
        }
        for analysis, subset in analyses.items():
            outcome = outcomes[analysis]
            subset = subset[subset[outcome].notna()].copy()
            if subset.empty or subset[outcome].nunique() < 2:
                continue
            for feature in loc_cols:
                present = subset[feature].fillna(0).astype(bool)
                positive = subset[outcome].astype(bool)
                a = int((present & positive).sum())
                b = int((present & ~positive).sum())
                c = int((~present & positive).sum())
                d = int((~present & ~positive).sum())
                odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="two-sided")
                rows.append({
                    "cancer_code": cancer,
                    "analysis": analysis,
                    "localization": feature,
                    "odds_ratio": odds,
                    "pvalue": pvalue,
                    "n_with_localization_positive": a,
                    "n_with_localization_negative": b,
                    "n_without_localization_positive": c,
                    "n_without_localization_negative": d,
                })
    results = pd.DataFrame(rows)
    if len(results):
        results["fdr"] = (
            results.groupby(["cancer_code", "analysis"], group_keys=False)["pvalue"]
            .transform(lambda x: multipletests(x, method="fdr_bh")[1])
        )
    return results


def main() -> None:
    master = pd.read_csv(DATA_PROCESSED / "master_table_augmented.csv")
    master["gene"] = master["gene"].astype(str).str.upper().str.strip()
    localization = pd.read_csv(DATA_PROCESSED / "hpa_subcellular_features_wide.csv")
    localization["gene"] = localization["gene"].astype(str).str.upper().str.strip()
    loc_cols = [column for column in localization if column.startswith("loc_")]

    cohort_tables = []
    cohort_summary = {}
    for directory, cancer_code in COHORTS.items():
        cohort_dir = ROOT / "data" / "cptac" / directory
        rna = pd.read_csv(cohort_dir / "matched_transcriptomics.csv.gz", index_col=0)
        protein = pd.read_csv(cohort_dir / "matched_proteomics.csv.gz", index_col=0)
        sample_info = pd.read_csv(cohort_dir / "matched_samples.csv")
        tumor_samples = pd.Index(
            sample_info.loc[
                sample_info["sample_kind_heuristic"] != "normal", "sample_id"
            ].astype(str)
        )
        tumor_samples = tumor_samples.intersection(rna.index).intersection(protein.index)
        metrics = gene_metrics(rna, protein, tumor_samples)
        metrics["cancer_code"] = cancer_code

        hpa = (
            master[master["cancer"].str.upper() == cancer_code]
            [["gene", "protein_detected", "detected_fraction", "total_patients"]]
            .drop_duplicates("gene")
        )
        metrics = metrics.merge(hpa, on="gene", how="inner", validate="one_to_one")
        metrics["hpa_negative"] = metrics["protein_detected"].eq(0)
        metrics["cptac_quantified"] = metrics["cptac_protein_coverage"].ge(0.80)
        metrics["cross_assay_disagreement"] = (
            metrics["hpa_negative"] & metrics["cptac_quantified"]
        )
        metrics["rna_protein_fdr"] = fdr(metrics["rna_protein_p"])
        cohort_tables.append(metrics)
        cohort_summary[cancer_code] = {
            "tumors": len(tumor_samples),
            "hpa_cptac_overlap_genes": len(metrics),
            "hpa_negative": int(metrics["hpa_negative"].sum()),
            "cptac_quantified": int(metrics["cptac_quantified"].sum()),
            "hpa_negative_cptac_quantified": int(metrics["cross_assay_disagreement"].sum()),
            "disagreement_rate_among_overlap": float(metrics["cross_assay_disagreement"].mean()),
            "disagreement_rate_among_hpa_negative": float(
                metrics.loc[metrics["hpa_negative"], "cptac_quantified"].mean()
            ),
            "median_rna_protein_spearman": float(metrics["rna_protein_spearman"].median()),
        }

    atlas = pd.concat(cohort_tables, ignore_index=True)
    atlas = atlas.merge(localization, on="gene", how="left", validate="many_to_one")
    atlas[loc_cols] = atlas[loc_cols].fillna(0).astype(int)

    recurrence = (
        atlas.groupby("gene")
        .agg(
            cancers_observed=("cancer_code", "nunique"),
            cancers_hpa_negative=("hpa_negative", "sum"),
            cancers_cptac_quantified=("cptac_quantified", "sum"),
            cancers_discordant=("cross_assay_disagreement", "sum"),
            median_cptac_coverage=("cptac_protein_coverage", "median"),
            median_hpa_detection_fraction=("detected_fraction", "median"),
            median_rna_protein_spearman=("rna_protein_spearman", "median"),
        )
        .reset_index()
        .sort_values(
            ["cancers_discordant", "median_cptac_coverage"], ascending=[False, False]
        )
    )
    enrichment = localization_tests(atlas, loc_cols)

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    atlas.to_csv(TABLES_DIR / "cross_cancer_hpa_cptac_atlas.csv", index=False)
    recurrence.to_csv(TABLES_DIR / "cross_cancer_disagreement_recurrence.csv", index=False)
    enrichment.to_csv(TABLES_DIR / "cross_cancer_localization_enrichment.csv", index=False)

    overview = {
        "cptac_quantified_threshold": 0.80,
        "cohorts": cohort_summary,
        "total_gene_cancer_pairs": len(atlas),
        "unique_genes": int(atlas["gene"].nunique()),
        "recurrent_disagreement": {
            str(n): int((recurrence["cancers_discordant"] >= n).sum())
            for n in range(1, len(COHORTS) + 1)
        },
        "localization_tests_fdr_significant": int((enrichment["fdr"] < 0.05).sum()),
    }
    (REPORTS_DIR / "cross_cancer_disagreement_overview.json").write_text(
        json.dumps(overview, indent=2), encoding="utf-8"
    )
    print(json.dumps(overview, indent=2))


if __name__ == "__main__":
    main()

