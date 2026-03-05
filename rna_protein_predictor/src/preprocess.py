# src/preprocess.py
import numpy as np
import pandas as pd
from .schema import Columns

C = Columns()

def compute_log2fc(tumor_tpm: pd.Series, normal_tpm: pd.Series, pseudocount: float = 1e-3) -> pd.Series:
    """
    log2((tumor + pc) / (normal + pc))
    """
    return np.log2((tumor_tpm + pseudocount) / (normal_tpm + pseudocount))

def protein_level_to_detected(level: pd.Series) -> pd.Series:
    """
    Map HPA categorical protein levels to binary detected target.
    detected=1 if High/Medium/Low; 0 if Not detected.
    """
    level_clean = level.astype(str).str.strip().str.lower()
    detected = ~level_clean.eq("not detected")
    return detected.astype(int)

def basic_clean(df: pd.DataFrame) -> pd.DataFrame:
    """
    Basic cleanup: standardize gene/cancer, drop obvious missing targets/features.
    """
    out = df.copy()
    out[C.gene] = out[C.gene].astype(str).str.strip().str.upper()
    out[C.cancer] = out[C.cancer].astype(str).str.strip().str.upper()
    return out