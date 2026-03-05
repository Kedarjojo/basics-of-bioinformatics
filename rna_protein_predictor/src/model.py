# src/model.py
import joblib
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from .schema import Columns, required_columns_for_v0

C = Columns()

def train_baseline_logreg(df: pd.DataFrame, random_state: int = 7):
    """
    Baseline: Protein_detected ~ RNA_log2FC
    """
    # validate minimal schema
    missing = [c for c in required_columns_for_v0() if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    X = df[[C.rna_log2fc]].copy()
    y = df[C.protein_detected].astype(int).copy()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=random_state
    )

    pipe = Pipeline(steps=[
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=200, class_weight="balanced"))
    ])

    pipe.fit(X_train, y_train)
    return pipe, (X_train, X_test, y_train, y_test)

def save_model(model, path):
    joblib.dump(model, path)

def load_model(path):
    return joblib.load(path)