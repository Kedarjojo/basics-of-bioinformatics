# src/analyze_discordance.py
import pandas as pd
from .config import DATA_PROCESSED
from .schema import Columns

C = Columns()

DATA = DATA_PROCESSED / "master_table_augmented.csv"

def main():
    df = pd.read_csv(DATA)

    discordant_rows = []

    for cancer, sub in df.groupby(C.cancer):
        threshold = sub["tumor_log2tpm_mean"].quantile(0.75)

        sub = sub.copy()
        sub["high_rna"] = sub["tumor_log2tpm_mean"] >= threshold
        sub["discordant"] = sub["high_rna"] & (sub[C.protein_detected] == 0)

        sub["cancer"] = cancer
        discordant_rows.append(sub)

    df_all = pd.concat(discordant_rows)

    print("Overall discordance rate among high-RNA genes:")
    print(
        df_all[df_all["high_rna"]]["discordant"].mean()
    )

    print("\nDiscordance by cancer:")
    print(
        df_all[df_all["high_rna"]]
        .groupby("cancer")["discordant"]
        .mean()
        .sort_values(ascending=False)
    )

    # Save discordant genes
    discordant_genes = df_all[df_all["discordant"]].copy()
    discordant_genes.to_csv(
        DATA_PROCESSED / "discordant_high_rna_genes.csv",
        index=False
    )

    print("\nTop discordant examples:")
    print(
        discordant_genes
        .sort_values("tumor_log2tpm_mean", ascending=False)
        .head(15)[["gene", "cancer", "tumor_log2tpm_mean", "rna_log2fc"]]
    )

    print("\nGene type distribution among discordant high-RNA genes:")
    print(
        discordant_genes["gene_type"]
        .value_counts(normalize=True)
        .head(10)
    )

    print("\nGene type distribution among concordant high-RNA genes:")
    concordant = df_all[(df_all["high_rna"]) & (~df_all["discordant"])]
    print(
        concordant["gene_type"]
        .value_counts(normalize=True)
        .head(10)
    )

    # Identify mitochondrial genes
    discordant_genes["mitochondrial"] = discordant_genes["gene"].str.startswith("MT-")
    concordant["mitochondrial"] = concordant["gene"].str.startswith("MT-")

    print("\nMitochondrial proportion:")
    print("Discordant:", discordant_genes["mitochondrial"].mean())
    print("Concordant:", concordant["mitochondrial"].mean())

    # Ribosomal genes
    discordant_genes["ribosomal"] = discordant_genes["gene"].str.startswith(("RPL", "RPS"))
    concordant["ribosomal"] = concordant["gene"].str.startswith(("RPL", "RPS"))

    print("\nRibosomal proportion:")
    print("Discordant:", discordant_genes["ribosomal"].mean())
    print("Concordant:", concordant["ribosomal"].mean())

    print("\nMean tumor expression:")
    print("Discordant:", discordant_genes["tumor_log2tpm_mean"].mean())
    print("Concordant:", concordant["tumor_log2tpm_mean"].mean())

    print("\nMedian tumor expression:")
    print("Discordant:", discordant_genes["tumor_log2tpm_mean"].median())
    print("Concordant:", concordant["tumor_log2tpm_mean"].median())

    print("\nMean gene length:")
    print("Discordant:", discordant_genes["gene_length"].mean())
    print("Concordant:", concordant["gene_length"].mean())

    gbm = df_all[(df_all["cancer"] == "GBM") & (df_all["high_rna"])]

    print("\nGBM discordance expression comparison:")
    print("Discordant mean:", gbm[gbm["discordant"]]["tumor_log2tpm_mean"].mean())
    print("Concordant mean:", gbm[~gbm["discordant"]]["tumor_log2tpm_mean"].mean())

    print("\nEffect sizes:")
    print("Expression difference:",
          concordant["tumor_log2tpm_mean"].mean() - discordant_genes["tumor_log2tpm_mean"].mean())

    print("Gene length difference:",
          concordant["gene_length"].mean() - discordant_genes["gene_length"].mean())

if __name__ == "__main__":
    main()