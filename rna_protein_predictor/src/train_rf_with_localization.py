# src/train_rf_with_localization.py
import pandas as pd
import joblib
from sklearn.model_selection import GroupShuffleSplit
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score

from .config import DATA_PROCESSED, MODELS_DIR
from .schema import Columns

C = Columns()

MASTER = DATA_PROCESSED / "master_table_augmented.csv"
LOC_WIDE = DATA_PROCESSED / "hpa_subcellular_features_wide.csv"
MODEL_OUT = MODELS_DIR / "rf_v2_with_localization.joblib"

BASE_FEATURES = [C.rna_log2fc, "tumor_log2tpm_mean", "gene_length"]

def main():
    df = pd.read_csv(MASTER)
    loc = pd.read_csv(LOC_WIDE)

    df["gene"] = df["gene"].astype(str).str.upper().str.strip()
    loc["gene"] = loc["gene"].astype(str).str.upper().str.strip()
    df = df.merge(loc, on="gene", how="left")

    loc_cols = [c for c in df.columns if c.startswith("loc_")]
    for c in loc_cols:
        df[c] = df[c].fillna(0).astype(int)

    FEATURES = BASE_FEATURES + loc_cols

    df = df.dropna(subset=BASE_FEATURES + [C.protein_detected])

    X = df[FEATURES]
    y = df[C.protein_detected].astype(int)
    groups = df[C.cancer]

    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=7)
    train_idx, test_idx = next(splitter.split(X, y, groups=groups))

    model = RandomForestClassifier(
        n_estimators=300,
        class_weight="balanced",
        random_state=7,
        n_jobs=-1,
    )
    model.fit(X.iloc[train_idx], y.iloc[train_idx])

    prob = model.predict_proba(X.iloc[test_idx])[:, 1]
    roc = roc_auc_score(y.iloc[test_idx], prob)
    pr = average_precision_score(y.iloc[test_idx], prob)

    print("Held-out cancers in TEST:", sorted(set(groups.iloc[test_idx])))
    print("ROC-AUC:", roc)
    print("PR-AUC:", pr)

    joblib.dump(model, MODEL_OUT)
    print("Saved:", MODEL_OUT)

    # importance summary
    imp = pd.Series(model.feature_importances_, index=FEATURES).sort_values(ascending=False)
    print("\nTop 20 importances:")
    print(imp.head(20).to_string())

if __name__ == "__main__":
    main()