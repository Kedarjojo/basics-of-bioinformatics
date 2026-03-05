# src/analyze_discordance_localization.py
import numpy as np
import pandas as pd
from .config import DATA_PROCESSED
from .schema import Columns

C = Columns()

MASTER = DATA_PROCESSED / "master_table_augmented.csv"
LOC_WIDE = DATA_PROCESSED / "hpa_subcellular_features_wide.csv"

OUT = DATA_PROCESSED / "discordance_localization_enrichment.csv"

def enrichment_table(df_high: pd.DataFrame, label_col: str, feature_cols: list[str]) -> pd.DataFrame:
    """
    Computes prevalence and enrichment for each binary feature:
      p(feature | discordant) / p(feature | concordant)
    """
    disc = df_high[df_high[label_col]].copy()
    conc = df_high[~df_high[label_col]].copy()

    rows = []
    for f in feature_cols:
        p_disc = disc[f].mean() if len(disc) else np.nan
        p_conc = conc[f].mean() if len(conc) else np.nan
        enr = (p_disc / p_conc) if (p_conc and p_conc > 0) else np.nan

        rows.append({
            "feature": f,
            "p_discordant": p_disc,
            "p_concordant": p_conc,
            "enrichment_ratio": enr,
            "n_discordant": len(disc),
            "n_concordant": len(conc),
        })

    out = pd.DataFrame(rows).sort_values("enrichment_ratio", ascending=False)
    return out

def main():
    df = pd.read_csv(MASTER)

    loc = pd.read_csv(LOC_WIDE)
    loc["gene"] = loc["gene"].astype(str).str.upper().str.strip()

    df["gene"] = df["gene"].astype(str).str.upper().str.strip()

    df = df.merge(loc, on="gene", how="left")

    # fill missing location features with 0
    loc_cols = [c for c in df.columns if c.startswith("loc_")]
    for c in loc_cols:
        df[c] = df[c].fillna(0).astype(int)

    # define high RNA + discordant within each cancer
    out_parts = []

    df_all = []
    for cancer, sub in df.groupby(C.cancer):
        sub = sub.copy()
        thr = sub["tumor_log2tpm_mean"].quantile(0.75)
        sub["high_rna"] = sub["tumor_log2tpm_mean"] >= thr
        sub["discordant"] = sub["high_rna"] & (sub[C.protein_detected].astype(int) == 0)
        df_all.append(sub)

        high = sub[sub["high_rna"]].copy()
        if len(high) == 0:
            continue

        enr = enrichment_table(high, "discordant", loc_cols)
        enr.insert(0, "cancer", cancer)
        out_parts.append(enr)

    df_all = pd.concat(df_all, ignore_index=True)
    enr_all = pd.concat(out_parts, ignore_index=True)

    # overall enrichment (pool cancers but still only high_rna)
    high_all = df_all[df_all["high_rna"]].copy()
    enr_overall = enrichment_table(high_all, "discordant", loc_cols)
    enr_overall.insert(0, "cancer", "ALL")

    final = pd.concat([enr_overall, enr_all], ignore_index=True)
    final.to_csv(OUT, index=False)

    print("Saved:", OUT)
    print("\nTop enriched localizations overall (discordant vs concordant, high-RNA only):")
    print(
        enr_overall
        .dropna(subset=["enrichment_ratio"])
        .sort_values("enrichment_ratio", ascending=False)
        .head(15)[["feature", "enrichment_ratio", "p_discordant", "p_concordant"]]
        .to_string(index=False)
    )

    print("\nTop enriched per cancer (top 5 each):")
    for cancer in sorted(df_all[C.cancer].unique()):
        sub = final[final["cancer"] == cancer].dropna(subset=["enrichment_ratio"])
        if len(sub) == 0:
            continue
        top5 = sub.sort_values("enrichment_ratio", ascending=False).head(5)
        print(f"\n{cancer}")
        print(top5[["feature", "enrichment_ratio", "p_discordant", "p_concordant"]].to_string(index=False))

if __name__ == "__main__":
    main()