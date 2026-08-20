"""Patient-level CPTAC discordance and HPA cross-assay triangulation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

from .config import DATA_PROCESSED, REPORTS_DIR, TABLES_DIR

ROOT = Path(__file__).resolve().parents[1]


def zscore(values: np.ndarray) -> np.ndarray:
    mean = np.nanmean(values)
    sd = np.nanstd(values, ddof=1)
    return (values - mean) / sd if np.isfinite(sd) and sd > 0 else np.full_like(values, np.nan)


def loo_residual(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Leave-one-out residuals for y ~ intercept + x using the hat diagonal."""
    design = np.column_stack([np.ones(len(x)), x])
    inverse = np.linalg.pinv(design.T @ design)
    coefficients = inverse @ design.T @ y
    fitted = design @ coefficients
    residual = y - fitted
    leverage = np.einsum("ij,jk,ik->i", design, inverse, design)
    denominator = np.clip(1 - leverage, 1e-8, None)
    return residual / denominator, float(coefficients[1]), float(coefficients[0])


def bh_adjust(pvalues: pd.Series) -> np.ndarray:
    valid = pvalues.notna()
    adjusted = np.full(len(pvalues), np.nan)
    if valid.any():
        adjusted[valid] = multipletests(pvalues[valid], method="fdr_bh")[1]
    return adjusted


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer", default="coad")
    parser.add_argument("--min-patients", type=int, default=30)
    parser.add_argument("--min-protein-coverage", type=float, default=0.50)
    parser.add_argument("--high-rna-z", type=float, default=1.0)
    parser.add_argument("--low-residual-z", type=float, default=-1.0)
    args = parser.parse_args()

    cohort_dir = ROOT / "data" / "cptac" / args.cancer
    rna = pd.read_csv(
        cohort_dir / "matched_transcriptomics.csv.gz", index_col=0
    )
    protein = pd.read_csv(
        cohort_dir / "matched_proteomics.csv.gz", index_col=0
    )
    shared_samples = rna.index.intersection(protein.index)
    shared_genes = rna.columns.intersection(protein.columns)
    rna = rna.loc[shared_samples, shared_genes]
    protein = protein.loc[shared_samples, shared_genes]

    rows: list[dict] = []
    patient_rows: list[dict] = []
    for gene in shared_genes:
        rna_values = pd.to_numeric(rna[gene], errors="coerce").to_numpy(float)
        protein_values = pd.to_numeric(protein[gene], errors="coerce").to_numpy(float)
        complete = np.isfinite(rna_values) & np.isfinite(protein_values)
        n = int(complete.sum())
        coverage = n / len(shared_samples)
        if n < args.min_patients or coverage < args.min_protein_coverage:
            continue

        x = zscore(rna_values[complete])
        y = zscore(protein_values[complete])
        if not np.isfinite(x).all() or not np.isfinite(y).all():
            continue
        rho, pvalue = spearmanr(x, y)
        residual, slope, intercept = loo_residual(x, y)
        residual_z = zscore(residual)
        high = x >= args.high_rna_z
        unexpected_low = high & (residual_z <= args.low_residual_z)
        low_rna = x <= -args.high_rna_z
        unexpected_high = low_rna & (residual_z >= -args.low_residual_z)

        rows.append({
            "gene": gene,
            "n_complete": n,
            "protein_coverage": coverage,
            "spearman_rho": float(rho),
            "spearman_p": float(pvalue),
            "rna_to_protein_slope": slope,
            "rna_to_protein_intercept": intercept,
            "n_high_rna": int(high.sum()),
            "n_high_rna_unexpected_low_protein": int(unexpected_low.sum()),
            "fraction_high_rna_unexpected_low_protein": (
                float(unexpected_low.sum() / high.sum()) if high.sum() else np.nan
            ),
            "median_residual_z_among_high_rna": (
                float(np.median(residual_z[high])) if high.sum() else np.nan
            ),
            "n_low_rna": int(low_rna.sum()),
            "n_low_rna_unexpected_high_protein": int(unexpected_high.sum()),
        })

        complete_samples = shared_samples[complete]
        for sample, rna_z, protein_z, resid_z in zip(complete_samples, x, y, residual_z):
            patient_rows.append({
                "sample_id": sample,
                "gene": gene,
                "rna_z": rna_z,
                "protein_z": protein_z,
                "loo_residual_z": resid_z,
                "high_rna_unexpected_low_protein": bool(
                    rna_z >= args.high_rna_z and resid_z <= args.low_residual_z
                ),
            })

    genes = pd.DataFrame(rows)
    genes["spearman_fdr"] = bh_adjust(genes["spearman_p"])

    master = pd.read_csv(DATA_PROCESSED / "master_table_augmented.csv")
    hpa = (
        master[master["cancer"].str.upper() == args.cancer.upper()]
        [["gene", "protein_detected", "detected_fraction", "total_patients"]]
        .drop_duplicates("gene")
    )
    hpa["gene"] = hpa["gene"].astype(str).str.upper().str.strip()
    genes["gene"] = genes["gene"].astype(str).str.upper().str.strip()
    genes = genes.merge(hpa, on="gene", how="left", validate="one_to_one")

    enough_high = genes["n_high_rna"] >= 5
    frequent_low = genes["fraction_high_rna_unexpected_low_protein"] >= 0.25
    well_observed = genes["protein_coverage"] >= 0.80
    hpa_negative = genes["protein_detected"] == 0
    hpa_positive = genes["protein_detected"] == 1
    tracks_rna = genes["spearman_rho"] >= 0.30

    genes["triangulation_class"] = "unclassified"
    genes.loc[hpa_negative & ~well_observed, "triangulation_class"] = "assay_ambiguous"
    genes.loc[
        hpa_negative & well_observed & enough_high & frequent_low,
        "triangulation_class",
    ] = "replicated_discordance_candidate"
    genes.loc[
        hpa_negative & well_observed & tracks_rna & ~frequent_low,
        "triangulation_class",
    ] = "possible_ihc_specific_nondetection"
    genes.loc[hpa_positive & well_observed, "triangulation_class"] = "cross_assay_detected"

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    genes = genes.sort_values(
        ["triangulation_class", "fraction_high_rna_unexpected_low_protein"],
        ascending=[True, False],
    )
    gene_path = TABLES_DIR / f"cptac_{args.cancer}_gene_discordance.csv"
    patient_path = TABLES_DIR / f"cptac_{args.cancer}_patient_discordance.csv.gz"
    genes.to_csv(gene_path, index=False)
    pd.DataFrame(patient_rows).to_csv(patient_path, index=False, compression="gzip")

    overview = {
        "cancer": args.cancer,
        "matched_patients": len(shared_samples),
        "shared_input_genes": len(shared_genes),
        "analyzed_genes": len(genes),
        "thresholds_are_provisional": True,
        "parameters": vars(args),
        "triangulation_counts": genes["triangulation_class"].value_counts().to_dict(),
        "median_spearman_rho": float(genes["spearman_rho"].median()),
        "hpa_overlap": int(genes["protein_detected"].notna().sum()),
    }
    report_path = REPORTS_DIR / f"cptac_{args.cancer}_discordance_overview.json"
    report_path.write_text(json.dumps(overview, indent=2), encoding="utf-8")
    print(json.dumps(overview, indent=2))
    print(f"Saved {gene_path}")
    print(f"Saved {patient_path}")


if __name__ == "__main__":
    main()

