# src/load_hpa.py
import pandas as pd
from pathlib import Path

from .config import DATA_RAW
from .schema import Columns
from .preprocess import basic_clean

C = Columns()

HPA_FILE = DATA_RAW / "hpa" / "cancer_data.tsv"

# Map HPA cancer names → TCGA codes
CANCER_MAP = {
    "breast cancer": "BRCA",
    "glioma": "GBM",
    "liver cancer": "LIHC",
    "pancreatic cancer": "PAAD",
    "thyroid cancer": "THCA",
    "prostate cancer": "PRAD",
    "colorectal cancer": "COAD",
}

def load_hpa_clean():
    df = pd.read_csv(HPA_FILE, sep="\t")

    # Rename columns to match schema
    df = df.rename(columns={
        "Gene name": C.gene,
        "Cancer": C.cancer,
    })

    # Keep only needed columns
    df = df[[C.gene, C.cancer, "High", "Medium", "Low", "Not detected"]]

    # Standardize
    df = basic_clean(df)
    df[C.cancer] = df[C.cancer].str.lower().map(CANCER_MAP)

    # Drop cancers not in mapping
    df = df.dropna(subset=[C.cancer])

    # Compute totals
    df["total_patients"] = (
        df["High"] + df["Medium"] + df["Low"] + df["Not detected"]
    )

    df["detected_count"] = df["High"] + df["Medium"] + df["Low"]

    # Majority rule
    df[C.protein_detected] = (
        df["detected_count"] > df["Not detected"]
    ).astype(int)

    # Optional: detection fraction
    df["detected_fraction"] = df["detected_count"] / df["total_patients"]

    return df[[C.gene, C.cancer, C.protein_detected, "detected_fraction", "total_patients"]]