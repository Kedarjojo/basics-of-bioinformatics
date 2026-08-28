#!/usr/bin/env python3
"""Download and audit independent reviewed-human UniProt annotations."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "outputs" / "tables" / "cptac_gene_predictability_atlas.csv"
ANNOTATION_DIR = ROOT / "data" / "annotations"
TABLE_DIR = ROOT / "outputs" / "tables"
REPORT_DIR = ROOT / "outputs" / "reports"

UNIPROT_STREAM_ENDPOINT = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_QUERY = "(proteome:UP000005640) AND (reviewed:true)"
UNIPROT_FIELDS = [
    "accession",
    "id",
    "gene_primary",
    "length",
    "cc_subcellular_location",
    "ft_signal",
    "ft_transmem",
    "protein_families",
]

EXPECTED_COLUMNS = [
    "accession",
    "entry_name",
    "gene",
    "protein_length",
    "subcellular_location_raw",
    "signal_peptide_raw",
    "transmembrane_raw",
    "protein_family_raw",
]

LOCATION_PATTERNS = {
    "uniprot_secreted": r"\bsecreted\b",
    "uniprot_extracellular": r"extracellular",
    "uniprot_plasma_membrane": r"cell membrane|plasma membrane",
    "uniprot_membrane_any": r"\bmembrane\b",
    "uniprot_cytosol": r"\bcytosol\b|\bcytoplasm\b",
    "uniprot_nucleus": r"\bnucleus\b|\bnuclear\b",
    "uniprot_mitochondrion": r"mitochond",
    "uniprot_endoplasmic_reticulum": r"endoplasmic reticulum",
    "uniprot_golgi": r"\bgolgi\b",
    "uniprot_lysosome": r"lysosom",
    "uniprot_peroxisome": r"peroxisom",
    "uniprot_cell_junction": r"cell junction|tight junction|adherens junction|desmosome",
}


def fetch_bytes(url: str, timeout: int) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "rna-protein-predictability/1.0"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def build_url() -> str:
    query = urllib.parse.urlencode(
        {
            "compressed": "true",
            "format": "tsv",
            "query": UNIPROT_QUERY,
            "fields": ",".join(UNIPROT_FIELDS),
        }
    )
    return f"{UNIPROT_STREAM_ENDPOINT}?{query}"


def download_uniprot(raw_path: Path, timeout: int, force: bool) -> tuple[bytes, str]:
    url = build_url()
    if raw_path.exists() and not force:
        return raw_path.read_bytes(), url
    payload = fetch_bytes(url, timeout)
    # UniProt returns gzip bytes when compressed=true.
    try:
        gzip.decompress(payload)
    except (gzip.BadGzipFile, EOFError) as error:
        preview = payload[:300].decode("utf-8", errors="replace")
        raise RuntimeError(f"UniProt response was not valid gzip: {preview}") from error
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_bytes(payload)
    return payload, url


def parse_uniprot(payload: bytes) -> pd.DataFrame:
    text = gzip.decompress(payload).decode("utf-8")
    frame = pd.read_csv(io.StringIO(text), sep="\t", dtype=str).fillna("")
    if frame.shape[1] != len(EXPECTED_COLUMNS):
        raise ValueError(
            "Unexpected UniProt schema. "
            f"Expected {len(EXPECTED_COLUMNS)} columns, received {frame.shape[1]}: "
            f"{frame.columns.tolist()}"
        )
    frame.columns = EXPECTED_COLUMNS
    frame["gene"] = frame["gene"].str.strip()
    frame["protein_length"] = pd.to_numeric(frame["protein_length"], errors="coerce")
    return frame


def count_transmembrane(value: str) -> int:
    return len(re.findall(r"\bTRANSMEM\b", value.upper()))


def collapse_gene_annotations(uniprot: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    usable = uniprot.loc[uniprot["gene"] != ""].copy()
    duplicate_counts = usable["gene"].value_counts()
    duplicate_genes = duplicate_counts.loc[duplicate_counts > 1].rename("uniprot_entries")
    duplicate_table = duplicate_genes.reset_index().rename(columns={"index": "gene"})

    usable["has_signal_peptide"] = usable["signal_peptide_raw"].str.strip().ne("")
    usable["transmembrane_domain_count"] = usable["transmembrane_raw"].map(
        count_transmembrane
    )
    location = usable["subcellular_location_raw"].str.lower()
    for column, pattern in LOCATION_PATTERNS.items():
        usable[column] = location.str.contains(pattern, regex=True, na=False)

    rows: list[dict] = []
    for gene, group in usable.groupby("gene", sort=True):
        row: dict[str, object] = {
            "gene": gene,
            "uniprot_annotation_available": True,
            "uniprot_entry_count": len(group),
            "uniprot_accessions": ";".join(sorted(group["accession"].unique())),
            "uniprot_reviewed": True,
            "protein_length": float(group["protein_length"].median()),
            "has_signal_peptide": bool(group["has_signal_peptide"].any()),
            "transmembrane_domain_count": int(
                group["transmembrane_domain_count"].max()
            ),
            "has_transmembrane_domain": bool(
                (group["transmembrane_domain_count"] > 0).any()
            ),
            "protein_family_available": bool(
                group["protein_family_raw"].str.strip().ne("").any()
            ),
            "protein_family_raw": " | ".join(
                sorted({x for x in group["protein_family_raw"] if x.strip()})
            ),
            "subcellular_location_raw": " | ".join(
                sorted({x for x in group["subcellular_location_raw"] if x.strip()})
            ),
        }
        for column in LOCATION_PATTERNS:
            row[column] = bool(group[column].any())
        rows.append(row)
    return pd.DataFrame(rows), duplicate_table


def join_atlas(atlas: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    merged = atlas[["gene"]].merge(
        annotations,
        on="gene",
        how="left",
        validate="one_to_one",
    )
    merged["uniprot_annotation_available"] = merged[
        "uniprot_annotation_available"
    ].map(lambda value: bool(value) if pd.notna(value) else False)
    boolean_columns = [
        "uniprot_reviewed",
        "has_signal_peptide",
        "has_transmembrane_domain",
        "protein_family_available",
        *LOCATION_PATTERNS.keys(),
    ]
    for column in boolean_columns:
        merged[column] = merged[column].map(
            lambda value: bool(value) if pd.notna(value) else False
        )
    merged["uniprot_entry_count"] = merged["uniprot_entry_count"].fillna(0).astype(int)
    merged["transmembrane_domain_count"] = (
        merged["transmembrane_domain_count"].fillna(0).astype(int)
    )
    for column in ("uniprot_accessions", "protein_family_raw", "subcellular_location_raw"):
        merged[column] = merged[column].fillna("")
    return merged


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, default=ATLAS)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if not args.atlas.exists():
        raise FileNotFoundError(f"Missing predictability atlas: {args.atlas}")
    atlas = pd.read_csv(args.atlas)
    if "gene" not in atlas.columns or atlas["gene"].duplicated().any():
        raise ValueError("Atlas must contain one unique row per gene")

    raw_path = ANNOTATION_DIR / "uniprot_reviewed_human.tsv.gz"
    payload, url = download_uniprot(raw_path, args.timeout, args.force)
    uniprot = parse_uniprot(payload)
    annotations, duplicate_table = collapse_gene_annotations(uniprot)
    atlas_annotations = join_atlas(atlas, annotations)

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    annotations.to_csv(TABLE_DIR / "uniprot_reviewed_human_gene_annotations.csv", index=False)
    atlas_annotations.to_csv(TABLE_DIR / "cptac_atlas_uniprot_annotations.csv", index=False)
    duplicate_table.to_csv(TABLE_DIR / "uniprot_duplicate_primary_gene_symbols.csv", index=False)

    available = atlas_annotations["uniprot_annotation_available"]
    report = {
        "analysis": "independent UniProt annotation coverage audit",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": "UniProtKB reviewed human reference proteome UP000005640",
        "query": UNIPROT_QUERY,
        "fields": UNIPROT_FIELDS,
        "request_url": url,
        "raw_sha256": hashlib.sha256(payload).hexdigest(),
        "raw_rows": len(uniprot),
        "raw_rows_with_primary_gene_symbol": int((uniprot["gene"] != "").sum()),
        "raw_unique_primary_gene_symbols": int(
            uniprot.loc[uniprot["gene"] != "", "gene"].nunique()
        ),
        "unique_annotated_genes": len(annotations),
        "duplicate_primary_gene_symbols": len(duplicate_table),
        "atlas_genes": len(atlas_annotations),
        "atlas_annotated_genes": int(available.sum()),
        "atlas_annotation_coverage": float(available.mean()),
        "coverage_gate_0_70_pass": bool(available.mean() >= 0.70),
        "feature_counts_among_atlas": {
            column: int(atlas_annotations[column].sum())
            for column in [
                "has_signal_peptide",
                "has_transmembrane_domain",
                *LOCATION_PATTERNS.keys(),
            ]
        },
        "missing_protein_length_among_annotated": int(
            atlas_annotations.loc[available, "protein_length"].isna().sum()
        ),
        "mapping_rule": (
            "Primary gene symbol; multiple reviewed entries collapsed by median length, "
            "maximum transmembrane count, and union of binary annotations"
        ),
    }
    (REPORT_DIR / "uniprot_annotation_coverage.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
