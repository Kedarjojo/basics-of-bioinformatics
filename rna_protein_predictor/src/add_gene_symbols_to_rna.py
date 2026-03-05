# src/add_gene_symbols_to_rna.py
import pandas as pd

from .config import DATA_RAW, DATA_PROCESSED
from .schema import Columns

C = Columns()

RNA_IN = DATA_PROCESSED / "rna_gene_cancer_table.csv"
RNA_OUT = DATA_PROCESSED / "rna_gene_cancer_table_with_symbols.csv"
HPA_CANCER_DATA = DATA_RAW / "hpa" / "cancer_data.tsv"

def main():
    rna = pd.read_csv(RNA_IN)

    # Strip version suffix: ENSG....2 -> ENSG....
    rna["ensembl_base"] = rna["ensembl_gene"].astype(str).str.split(".").str[0]

    hpa = pd.read_csv(HPA_CANCER_DATA, sep="\t")
    mapping = (
        hpa[["Gene", "Gene name"]]
        .drop_duplicates()
        .rename(columns={"Gene": "ensembl_base", "Gene name": C.gene})
    )

    mapping["ensembl_base"] = mapping["ensembl_base"].astype(str).str.split(".").str[0]
    mapping[C.gene] = mapping[C.gene].astype(str).str.upper().str.strip()

    merged = rna.merge(mapping, on="ensembl_base", how="left")

    # Drop rows without symbols (will be relatively few)
    before = len(merged)
    merged = merged.dropna(subset=[C.gene])
    after = len(merged)

    # Keep tidy columns
    out = merged[[C.gene, C.cancer, "ensembl_gene", "tumor_log2tpm_mean", "normal_log2tpm_mean", C.rna_log2fc]].copy()

    out.to_csv(RNA_OUT, index=False)

    print("Saved:", RNA_OUT)
    print(f"Rows kept: {after:,} / {before:,} (dropped {before-after:,} without symbols)")
    print("Unique genes:", out[C.gene].nunique())
    print("Cancers:", sorted(out[C.cancer].unique()))

if __name__ == "__main__":
    main()