# src/train_rf.py
import pandas as pd
import joblib

from sklearn.model_selection import GroupShuffleSplit
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, classification_report

from .config import DATA_PROCESSED, MODELS_DIR
from .schema import Columns

C = Columns()

DATA = DATA_PROCESSED / "master_table_augmented.csv"
MODEL_OUT = MODELS_DIR / "rf_v1.joblib"

FEATURES = [C.rna_log2fc, "tumor_log2tpm_mean", "gene_length", "protein_coding"]

def main():
    df = pd.read_csv(DATA)
    df = df.dropna(subset=FEATURES + [C.protein_detected])

    X = df[FEATURES]
    y = df[C.protein_detected].astype(int)
    groups = df[C.cancer]

    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=7)
    train_idx, test_idx = next(splitter.split(X, y, groups=groups))

    model = RandomForestClassifier(
        n_estimators=200,
        max_depth=None,
        class_weight="balanced",
        random_state=7,
        n_jobs=-1,
    )

    model.fit(X.iloc[train_idx], y.iloc[train_idx])

    prob = model.predict_proba(X.iloc[test_idx])[:, 1]
    pred = model.predict(X.iloc[test_idx])

    print("Held-out cancers in TEST:", sorted(set(groups.iloc[test_idx])))
    print("ROC-AUC:", roc_auc_score(y.iloc[test_idx], prob))
    print("PR-AUC:", average_precision_score(y.iloc[test_idx], prob))
    print(classification_report(y.iloc[test_idx], pred, digits=3))

    joblib.dump(model, MODEL_OUT)
    print("Saved:", MODEL_OUT)

    # Feature importance
    importances = pd.Series(model.feature_importances_, index=FEATURES)
    print("\nFeature Importances:")
    print(importances.sort_values(ascending=False))

if __name__ == "__main__":
    main()