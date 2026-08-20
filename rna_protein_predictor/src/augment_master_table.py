# src/augment_master_table.py
import pandas as pd

from .config import DATA_PROCESSED
from .schema import Columns

C = Columns()

MASTER_IN = DATA_PROCESSED / "master_table.csv"
GENE_FEAT = DATA_PROCESSED / "gene_basic_features.csv"
MASTER_OUT = DATA_PROCESSED / "master_table_augmented.csv"

def main():
    df = pd.read_csv(MASTER_IN)
    gf = pd.read_csv(GENE_FEAT)

    # protein_coding flag
    gf["protein_coding"] = (gf["gene_type"] == "protein_coding").astype(int)

    # Ensembl identifiers are stable join keys; symbols are display labels.
    df[C.gene] = df[C.gene].astype(str).str.upper().str.strip()
    gf["gene"] = gf["gene"].astype(str).str.upper().str.strip()
    df["ensembl_base"] = df["ensembl_base"].astype(str).str.split(".").str[0]
    gf["ensembl_base"] = gf["ensembl_base"].astype(str).str.split(".").str[0]

    merged = df.merge(
        gf[["ensembl_base", "gene_length", "gene_type", "protein_coding"]],
        on="ensembl_base",
        how="left",
        validate="many_to_one",
    )

    print("Before:", df.shape, "After:", merged.shape)
    print("Missing gene_length:", merged["gene_length"].isna().mean())

    merged.to_csv(MASTER_OUT, index=False)
    print("Saved:", MASTER_OUT)

if __name__ == "__main__":
    main()
