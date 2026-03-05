import pandas as pd

from src.schema import Columns
from src.preprocess import protein_level_to_detected, basic_clean
from src.model import train_baseline_logreg
from src.evaluate import evaluate_binary
from src.config import MODELS_DIR

C = Columns()

# Minimal dummy dataset (replace with real merged data later)
df = pd.DataFrame({
    C.gene: ["EMC9"]*6 + ["DBNDD2"]*6 + ["TP53"]*6 + ["GAPDH"]*6,
    C.cancer: ["DLBC","THYM","BRCA","LUAD","PAAD","COAD"]*4,
    C.rna_log2fc: [1.6,1.2,0.3,0.2,0.7,0.4,  0.9,1.1,0.8,0.6,0.5,0.7,
                  0.2,0.1,0.2,0.0,0.1,0.3, -0.1,0.0,0.1,-0.2,0.0,-0.1],
    C.hpa_protein_level: (
        ["Not detected","Not detected","Low","Low","Not detected","Low"] +   # EMC9
        ["High","Medium","High","Low","Medium","High"] +                    # DBNDD2
        ["High","High","High","High","High","High"] +                       # TP53
        ["High","High","High","High","High","High"]                         # GAPDH
    )
})

df[C.protein_detected] = protein_level_to_detected(df[C.hpa_protein_level])
df = basic_clean(df)

model, (X_train, X_test, y_train, y_test) = train_baseline_logreg(df)
metrics = evaluate_binary(model, X_test, y_test)

print("ROC-AUC:", metrics["roc_auc"])
print("PR-AUC:", metrics["pr_auc"])
print(metrics["report"])

# save model
model_path = MODELS_DIR / "baseline_logreg.joblib"
import joblib
joblib.dump(model, model_path)
print("Saved:", model_path)