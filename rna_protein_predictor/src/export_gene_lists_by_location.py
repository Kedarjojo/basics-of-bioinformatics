# src/export_gene_lists_by_location.py
import pandas as pd
from .config import DATA_PROCESSED

MASTER = DATA_PROCESSED / "master_table_augmented.csv"
LOC_WIDE = DATA_PROCESSED / "hpa_subcellular_features_wide.csv"
OUT = DATA_PROCESSED / "discordant_gene_lists_by_location.json"

def main():
    df = pd.read_csv(MASTER)
    loc = pd.read_csv(LOC_WIDE)

    df["gene"] = df["gene"].astype(str).str.upper().str.strip()
    loc["gene"] = loc["gene"].astype(str).str.upper().str.strip()

    df = df.merge(loc, on="gene", how="left")
    loc_cols = [c for c in df.columns if c.startswith("loc_")]
    for c in loc_cols:
        df[c] = df[c].fillna(0).astype(int)

    # define high RNA + discordant within each cancer
    out = {}

    for cancer, sub in df.groupby("cancer"):
        sub = sub.copy()
        thr = sub["tumor_log2tpm_mean"].quantile(0.75)
        sub["high_rna"] = sub["tumor_log2tpm_mean"] >= thr
        sub["discordant"] = sub["high_rna"] & (sub["protein_detected"].astype(int) == 0)

        disc = sub[sub["discordant"]].copy()

        out[cancer] = {}
        for c in loc_cols:
            genes = (
                disc.loc[disc[c] == 1, "gene"]
                .dropna()
                .unique()
                .tolist()
            )
            if len(genes) >= 15:  # avoid tiny lists
                out[cancer][c.replace("loc_", "")] = sorted(genes)

    import json
    OUT.write_text(json.dumps(out, indent=2))
    print("Saved:", OUT)

if __name__ == "__main__":
    main()