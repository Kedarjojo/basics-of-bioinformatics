# src/analyze_discordance_localization_stats.py
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from .config import DATA_PROCESSED
from .schema import Columns

C = Columns()
MASTER = DATA_PROCESSED / "master_table_augmented.csv"
LOC_WIDE = DATA_PROCESSED / "hpa_subcellular_features_wide.csv"
OUT = DATA_PROCESSED / "discordance_localization_stats.csv"

def fisher_for_feature(df_high, label_col, feature):
    # contingency: [[a,b],[c,d]] where a = discordant with feature, b = discordant without feature
    disc = df_high[df_high[label_col]]
    conc = df_high[~df_high[label_col]]
    a = int(disc[feature].sum())
    b = int(len(disc) - a)
    c = int(conc[feature].sum())
    d = int(len(conc) - c)
    # need to handle degenerate cases
    table = np.array([[a, b],[c,d]])
    try:
        oddsratio, p = fisher_exact(table, alternative="two-sided")
    except Exception:
        oddsratio, p = np.nan, 1.0
    return a, b, c, d, oddsratio, p

def main():
    df = pd.read_csv(MASTER)
    loc = pd.read_csv(LOC_WIDE)
    loc["gene"] = loc["gene"].astype(str).str.upper().str.strip()
    df["gene"] = df["gene"].astype(str).str.upper().str.strip()
    df = df.merge(loc, on="gene", how="left")

    loc_cols = [c for c in df.columns if c.startswith("loc_")]
    for c in loc_cols:
        df[c] = df[c].fillna(0).astype(int)

    rows = []
    df_all = []
    for cancer, sub in df.groupby(C.cancer):
        sub = sub.copy()
        thr = sub["tumor_log2tpm_mean"].quantile(0.75)
        sub["high_rna"] = sub["tumor_log2tpm_mean"] >= thr
        sub["discordant"] = sub["high_rna"] & (sub[C.protein_detected].astype(int) == 0)

        high = sub[sub["high_rna"]].copy()
        if len(high) == 0:
            continue

        for f in loc_cols:
            a,b,c,d,oratio,p = fisher_for_feature(high, "discordant", f)
            rows.append({
                "cancer": cancer,
                "feature": f,
                "n_disc_with": a,
                "n_disc_total": len(high[high["discordant"]]),
                "n_conc_with": c,
                "n_conc_total": len(high[~high["discordant"]]),
                "oddsratio": oratio,
                "pvalue": p
            })

    df_stats = pd.DataFrame(rows)
    # FDR within each cancer separately
    out_frames = []
    for cancer, sub in df_stats.groupby("cancer"):
        sub = sub.copy()
        rej, p_adj, _, _ = multipletests(sub["pvalue"].fillna(1.0), method="fdr_bh")
        sub["p_adj"] = p_adj
        sub["significant"] = rej
        out_frames.append(sub)
    final = pd.concat(out_frames, ignore_index=True)
    final = final.sort_values(["cancer", "p_adj"])
    final.to_csv(OUT, index=False)
    print("Saved:", OUT)
    print("Top significant (per cancer):")
    for cancer, sub in final.groupby("cancer"):
        top = sub[sub["p_adj"] < 0.05].sort_values("oddsratio", ascending=False).head(5)
        if len(top):
            print("\n", cancer)
            print(top[["feature","oddsratio","p_adj","n_disc_with","n_disc_total","n_conc_with","n_conc_total"]].to_string(index=False))

if __name__ == "__main__":
    main()