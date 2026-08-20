#!/usr/bin/env python3
"""Download and audit matched CPTAC RNA/protein data for one pilot cancer."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import cptac
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "cptac"
REPORTS = ROOT / "outputs" / "reports"

DATASETS = {
    "brca": ("Brca", "BRCA"),
    "coad": ("Coad", "COAD"),
    "gbm": ("Gbm", "GBM"),
    "pdac": ("Pdac", "PAAD"),
}


def patient_index(frame: pd.DataFrame) -> pd.Index:
    """Return CPTAC patient/sample labels without discarding tumor/normal suffixes."""
    if isinstance(frame.index, pd.MultiIndex):
        preferred = [
            name for name in frame.index.names
            if name and name.lower() in {"patient_id", "sample_id"}
        ]
        values = frame.index.get_level_values(preferred[0] if preferred else 0)
    else:
        values = frame.index
    return pd.Index(values.astype(str), name="sample_id")


def flatten_genes(frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse transcript/protein identifiers to one median column per gene symbol."""
    copy = frame.copy()
    if isinstance(copy.columns, pd.MultiIndex):
        genes = copy.columns.get_level_values(0).astype(str)
    else:
        genes = copy.columns.astype(str)
    copy.columns = genes
    copy.index = patient_index(copy)
    copy = copy.apply(pd.to_numeric, errors="coerce")
    # Transcripts and protein isoforms may map to the same symbol.
    copy = copy.T.groupby(level=0, sort=True).median().T
    return copy


def sample_kind(sample: str) -> str:
    upper = sample.upper()
    if upper.endswith(".N") or upper.endswith("-N") or upper.endswith("_N"):
        return "normal"
    return "tumor_or_unspecified"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer", choices=sorted(DATASETS), default="coad")
    parser.add_argument("--rna-source", default="washu")
    parser.add_argument("--protein-source", default="umich")
    args = parser.parse_args()

    class_name, cancer_code = DATASETS[args.cancer]
    dataset_class = getattr(cptac, class_name)
    print(f"Loading CPTAC {args.cancer}; the first call may download data...")
    dataset = dataset_class()
    print(dataset.list_data_sources())

    rna_raw = dataset.get_transcriptomics(args.rna_source)
    protein_raw = dataset.get_proteomics(args.protein_source)
    rna = flatten_genes(rna_raw)
    protein = flatten_genes(protein_raw)

    shared_samples = rna.index.intersection(protein.index)
    shared_genes = rna.columns.intersection(protein.columns)
    rna = rna.loc[shared_samples, shared_genes]
    protein = protein.loc[shared_samples, shared_genes]

    cancer_dir = OUT / args.cancer
    cancer_dir.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    rna.to_csv(cancer_dir / "matched_transcriptomics.csv.gz", compression="gzip")
    protein.to_csv(cancer_dir / "matched_proteomics.csv.gz", compression="gzip")

    samples = pd.DataFrame({"sample_id": shared_samples.astype(str)})
    samples["sample_kind_heuristic"] = samples["sample_id"].map(sample_kind)
    samples.to_csv(cancer_dir / "matched_samples.csv", index=False)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "cptac_package_version": getattr(cptac, "__version__", "unknown"),
        "cancer": args.cancer,
        "cancer_code": cancer_code,
        "rna_source": args.rna_source,
        "protein_source": args.protein_source,
        "rna_raw_shape": list(rna_raw.shape),
        "protein_raw_shape": list(protein_raw.shape),
        "matched_shape": [len(shared_samples), len(shared_genes)],
        "sample_kind_counts": samples["sample_kind_heuristic"].value_counts().to_dict(),
        "sample_examples": samples["sample_id"].head(12).tolist(),
        "rna_missing_fraction": float(rna.isna().mean().mean()),
        "protein_missing_fraction": float(protein.isna().mean().mean()),
    }
    manifest_path = REPORTS / f"cptac_{args.cancer}_pilot_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    print(f"Saved matched matrices under {cancer_dir.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

