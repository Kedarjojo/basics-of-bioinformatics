# src/build_master_table.py
import pandas as pd
from pathlib import Path

from .config import DATA_PROCESSED
from .schema import Columns
from .load_hpa import load_hpa_clean

C = Columns()

RNA_PATH = DATA_PROCESSED / "rna_gene_cancer_table_with_symbols.csv"
OUT_PATH = DATA_PROCESSED / "master_table.csv"

def build_master_table() -> pd.DataFrame:
    # HPA is already gene symbol + cancer code + protein_detected
    hpa = load_hpa_clean()

    # RNA table currently has ensembl_gene and cancer + rna_log2fc
    rna = pd.read_csv(RNA_PATH)

    # Standardize
    hpa[C.gene] = hpa[C.gene].astype(str).str.upper().str.strip()
    hpa[C.cancer] = hpa[C.cancer].astype(str).str.upper().str.strip()
    rna[C.gene] = rna[C.gene].astype(str).str.upper().str.strip()
    rna[C.cancer] = rna[C.cancer].astype(str).str.upper().str.strip()

    merged = pd.merge(
        hpa,
        rna[[C.gene, C.cancer, C.rna_log2fc, "tumor_log2tpm_mean", "normal_log2tpm_mean"]],
        on=[C.gene, C.cancer],
        how="inner",
    )

    # Drop any weird missing values
    merged = merged.dropna(subset=[C.rna_log2fc, C.protein_detected])

    return merged

def main():
    df = build_master_table()
    df.to_csv(OUT_PATH, index=False)
    print("Saved:", OUT_PATH)
    print("Rows:", len(df))
    print("Genes:", df[C.gene].nunique(), "Cancers:", df[C.cancer].nunique())
    print("Target counts:\n", df[C.protein_detected].value_counts())

if __name__ == "__main__":
    main()