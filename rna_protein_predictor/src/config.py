# src/config.py
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
OUTPUTS = PROJECT_ROOT / "outputs"
MODELS_DIR = OUTPUTS / "models"
FIGURES_DIR = OUTPUTS / "figures"
TABLES_DIR = OUTPUTS / "tables"
REPORTS_DIR = OUTPUTS / "reports"

for p in [
    DATA_RAW,
    DATA_PROCESSED,
    OUTPUTS,
    MODELS_DIR,
    FIGURES_DIR,
    TABLES_DIR,
    REPORTS_DIR,
]:
    p.mkdir(parents=True, exist_ok=True)
