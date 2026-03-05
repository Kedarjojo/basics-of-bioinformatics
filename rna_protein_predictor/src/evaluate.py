# src/evaluate.py
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, classification_report

def evaluate_binary(model, X_test: pd.DataFrame, y_test: pd.Series) -> dict:
    y_test = pd.Series(y_test).astype(int)
    prob = model.predict_proba(X_test)[:, 1]
    pred = model.predict(X_test)

    out = {
        "n_test": int(len(y_test)),
        "class_counts": y_test.value_counts().to_dict(),
        "pr_auc": average_precision_score(y_test, prob),
        "report": classification_report(y_test, pred, digits=3, zero_division=0),
    }

    # ROC-AUC only defined if both classes present
    if y_test.nunique() == 2:
        out["roc_auc"] = roc_auc_score(y_test, prob)
    else:
        out["roc_auc"] = float("nan")
        out["roc_auc_note"] = "ROC-AUC undefined because y_test contains only one class."

    return out