# src/make_localization_report_tables.py
import pandas as pd
from .config import DATA_PROCESSED

STATS = DATA_PROCESSED / "discordance_localization_stats.csv"
ENRICH = DATA_PROCESSED / "discordance_localization_enrichment.csv"

OUT_TOP = DATA_PROCESSED / "discordance_localization_top.csv"
OUT_WEB = DATA_PROCESSED / "discordance_localization_top.json"

def clean_feature_name(f: str) -> str:
    return f.replace("loc_", "").replace("_", " ")

def main():
    stats = pd.read_csv(STATS)
    enrich = pd.read_csv(ENRICH)

    # Merge enrichment ratio into stats (same cancer+feature)
    merged = stats.merge(
        enrich[["cancer", "feature", "enrichment_ratio", "p_discordant", "p_concordant"]],
        on=["cancer", "feature"],
        how="left",
    )

    merged["location"] = merged["feature"].map(clean_feature_name)

    # Choose “top” as: significant + large OR away from 1
    sig = merged[merged["p_adj"] < 0.05].copy()

    # Rank by |log(OR)|
    sig["abs_log_or"] = (sig["oddsratio"].astype(float)).apply(lambda x: abs(0 if x <= 0 else __import__("math").log(x)))
    sig = sig.sort_values(["cancer", "abs_log_or"], ascending=[True, False])

    # Take top 8 per cancer
    top = sig.groupby("cancer").head(8).copy()

    # Keep nice columns
    top = top[[
        "cancer", "location",
        "oddsratio", "p_adj",
        "enrichment_ratio", "p_discordant", "p_concordant",
        "n_disc_with", "n_disc_total", "n_conc_with", "n_conc_total"
    ]]

    top.to_csv(OUT_TOP, index=False)

    # Also export a JSON for webpage use
    top.to_json(OUT_WEB, orient="records", indent=2)

    print("Saved:", OUT_TOP)
    print("Saved:", OUT_WEB)
    print("\nPreview:")
    print(top.head(15).to_string(index=False))

if __name__ == "__main__":
    main()