# src/find_unexpected_missing_proteins.py
import pandas as pd
import joblib

from sklearn.model_selection import GroupShuffleSplit

from .config import DATA_PROCESSED, MODELS_DIR
from .schema import Columns

C = Columns()

MASTER = DATA_PROCESSED / "master_table_augmented.csv"
LOC_WIDE = DATA_PROCESSED / "hpa_subcellular_features_wide.csv"
MODEL_PATH = MODELS_DIR / "rf_v2_with_localization.joblib"  # your latest model

OUT_ALL = DATA_PROCESSED / "unexpected_missing_proteins_top.csv"
OUT_PER_CANCER = DATA_PROCESSED / "unexpected_missing_proteins_top_per_cancer.csv"

BASE_FEATURES = [C.rna_log2fc, "tumor_log2tpm_mean", "gene_length"]

def main():
    df = pd.read_csv(MASTER)
    loc = pd.read_csv(LOC_WIDE)

    # normalize gene symbols
    df["gene"] = df["gene"].astype(str).str.upper().str.strip()
    loc["gene"] = loc["gene"].astype(str).str.upper().str.strip()

    # merge localization features
    df = df.merge(loc, on="gene", how="left")
    loc_cols = [c for c in df.columns if c.startswith("loc_")]
    for c in loc_cols:
        df[c] = df[c].fillna(0).astype(int)

    FEATURES = BASE_FEATURES + loc_cols

    # drop missing core features/target
    df = df.dropna(subset=BASE_FEATURES + [C.protein_detected]).copy()
    df[C.protein_detected] = df[C.protein_detected].astype(int)

    # define high-RNA within each cancer (top quartile)
    df["high_rna"] = False
    for cancer, sub_idx in df.groupby(C.cancer).groups.items():
        sub = df.loc[sub_idx]
        thr = sub["tumor_log2tpm_mean"].quantile(0.75)
        df.loc[sub_idx, "high_rna"] = df.loc[sub_idx, "tumor_log2tpm_mean"] >= thr

    # consistent split: hold out cancers like you did in training/eval
    X = df[FEATURES]
    y = df[C.protein_detected]
    groups = df[C.cancer].astype(str)

    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=7)
    train_idx, test_idx = next(splitter.split(X, y, groups=groups))

    test = df.iloc[test_idx].copy()

    model = joblib.load(MODEL_PATH)
    test["pred_prob_detected"] = model.predict_proba(test[FEATURES])[:, 1]

    # "unexpected missing": high RNA, predicted detected, actually not detected
    # you can tune prob threshold; 0.80 is a good start
    PROB_THR = 0.80
    unexpected = test[
        (test["high_rna"]) &
        (test["pred_prob_detected"] >= PROB_THR) &
        (test[C.protein_detected] == 0)
    ].copy()

    print("Held-out cancers in TEST:", sorted(set(test[C.cancer])))
    print("Unexpected missing (count):", len(unexpected))
    print("Rate among high-RNA test rows:",
          len(unexpected) / max(1, test[test["high_rna"]].shape[0]))

    # sort by confidence then RNA abundance
    unexpected = unexpected.sort_values(
        ["pred_prob_detected", "tumor_log2tpm_mean"],
        ascending=[False, False]
    )

    # pick useful columns for reporting (make ensembl column optional)
    possible_ens = ["ensembl_gene", "ensembl", "Gene", "ensembl_id", "ensembl_base"]
    ens_col = next((c for c in possible_ens if c in test.columns), None)

    base_cols = [
        "gene", C.cancer,
        "tumor_log2tpm_mean", "normal_log2tpm_mean", C.rna_log2fc,
        C.protein_detected, "pred_prob_detected"
    ]

    if ens_col is not None and ens_col not in base_cols:
        base_cols.insert(2, ens_col)  # put it after cancer

    cols = base_cols + [c for c in loc_cols if c in unexpected.columns]

    # overall top 200
    unexpected[cols].head(200).to_csv(OUT_ALL, index=False)
    print("Saved:", OUT_ALL)

    # top 50 per cancer
    top_per = (
        unexpected[cols]
        .groupby(C.cancer, group_keys=False)
        .head(50)
    )
    top_per.to_csv(OUT_PER_CANCER, index=False)
    print("Saved:", OUT_PER_CANCER)

    # quick summary by cancer
    print("\nCounts by cancer:")
    print(unexpected.groupby(C.cancer).size().sort_values(ascending=False))

    print("\nTop 15 overall:")
    print(
        unexpected[["gene", C.cancer, "tumor_log2tpm_mean", C.rna_log2fc, "pred_prob_detected"]]
        .head(15)
        .to_string(index=False)
    )

if __name__ == "__main__":
    main()