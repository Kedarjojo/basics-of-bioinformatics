# src/train_baseline.py
import pandas as pd

from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, classification_report

from .config import DATA_PROCESSED, MODELS_DIR
from .schema import Columns
import joblib

C = Columns()

MASTER = DATA_PROCESSED / "master_table.csv"
MODEL_OUT = MODELS_DIR / "baseline_logreg_rna_log2fc.joblib"

def main():
    df = pd.read_csv(MASTER)

    # Features / target
    X = df[[C.rna_log2fc, "tumor_log2tpm_mean"]].copy()
    y = df[C.protein_detected].astype(int).copy()

    # IMPORTANT: evaluate generalization across cancers
    groups = df[C.cancer].astype(str)

    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=7)
    train_idx, test_idx = next(splitter.split(X, y, groups=groups))

    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    model = Pipeline(steps=[
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=500, class_weight="balanced"))
    ])

    model.fit(X_train, y_train)

    prob = model.predict_proba(X_test)[:, 1]
    pred = model.predict(X_test)

    roc = roc_auc_score(y_test, prob)
    pr = average_precision_score(y_test, prob)

    print("Held-out cancers in TEST:", sorted(set(groups.iloc[test_idx])))
    print("ROC-AUC:", roc)
    print("PR-AUC:", pr)
    print(classification_report(y_test, pred, digits=3))

    joblib.dump(model, MODEL_OUT)
    print("Saved model:", MODEL_OUT)

if __name__ == "__main__":
    main()