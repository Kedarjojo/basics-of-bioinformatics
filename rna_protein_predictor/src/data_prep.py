"""Shared, non-destructive preparation for audited model evaluation."""

from __future__ import annotations

import pandas as pd

from .config import DATA_PROCESSED

MASTER = DATA_PROCESSED / "master_table_augmented.csv"
LOCALIZATION = DATA_PROCESSED / "hpa_subcellular_features_wide.csv"
KEY = ["ensembl_base", "cancer"]


def load_model_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return a deterministic provisional table and its duplicate diagnostics.

    Duplicate resolution is provisional because the definitive fix is an
    Ensembl-ID join. For model evaluation, numeric conflicts are represented by
    their median and categorical conflicts by the first sorted non-null value.
    The raw source is never modified.
    """
    master = pd.read_csv(MASTER)
    loc = pd.read_csv(LOCALIZATION)
    for frame in (master, loc):
        frame["gene"] = frame["gene"].astype(str).str.upper().str.strip()
    master["cancer"] = master["cancer"].astype(str).str.upper().str.strip()

    duplicated = master[master.duplicated(KEY, keep=False)].copy()
    numeric = [c for c in master.select_dtypes(include="number") if c not in KEY]
    categorical = [c for c in master.columns if c not in KEY + numeric]
    aggregations = {c: "median" for c in numeric}
    aggregations.update({c: _first_sorted for c in categorical})
    clean = master.groupby(KEY, as_index=False, sort=True).agg(aggregations)

    clean = clean.merge(loc, on="gene", how="left", validate="many_to_one")
    loc_cols = [c for c in clean if c.startswith("loc_")]
    clean["localization_available"] = clean[loc_cols].notna().any(axis=1).astype(int)
    clean[loc_cols] = clean[loc_cols].fillna(0).astype(int)
    return clean, duplicated


def _first_sorted(values: pd.Series):
    non_null = sorted({str(v) for v in values.dropna()})
    return non_null[0] if non_null else pd.NA

