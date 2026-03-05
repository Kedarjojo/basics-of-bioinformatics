# src/load_rna.py
from pathlib import Path
import pandas as pd
import numpy as np

from .config import DATA_RAW
from .schema import Columns

C = Columns()

PHENO_FILE = DATA_RAW / "rnaseq" / "TcgaTargetGTEX_phenotype.txt"
# Change this to your actual gene TPM filename:
EXPR_FILE = DATA_RAW / "rnaseq" / "TcgaTargetGtex_rsem_gene_tpm.gz"

TCGA_DISEASE_TO_CODE = {
    "Breast Invasive Carcinoma": "BRCA",
    "Colon Adenocarcinoma": "COAD",
    "Glioblastoma Multiforme": "GBM",
    "Liver Hepatocellular Carcinoma": "LIHC",
    "Pancreatic Adenocarcinoma": "PAAD",
    "Prostate Adenocarcinoma": "PRAD",
    "Thyroid Carcinoma": "THCA",
}

CANCER_TO_GTEX_SITE = {
    "BRCA": "Breast",
    "COAD": "Colon",
    "GBM": "Brain",
    "LIHC": "Liver",
    "PAAD": "Pancreas",
    "PRAD": "Prostate",
    "THCA": "Thyroid",
}

def load_pheno() -> pd.DataFrame:
    pheno = pd.read_csv(PHENO_FILE, sep="\t", encoding="latin-1")
    return pheno

def select_samples(pheno: pd.DataFrame):
    # TCGA tumors for selected cancers
    tcga = pheno[(pheno["_study"] == "TCGA") & (pheno["_sample_type"] == "Primary Tumor")].copy()
    tcga = tcga[tcga["primary disease or tissue"].isin(TCGA_DISEASE_TO_CODE.keys())].copy()
    tcga["cancer_code"] = tcga["primary disease or tissue"].map(TCGA_DISEASE_TO_CODE)

    # GTEx normals for matching sites
    gtex = pheno[(pheno["_study"] == "GTEX") & (pheno["_sample_type"] == "Normal Tissue")].copy()

    # Build per-cancer GTEx sample lists
    gtex_samples_by_cancer = {}
    for code, site in CANCER_TO_GTEX_SITE.items():
        if site.lower() == "brain":
            mask = gtex["_primary_site"].astype(str).str.contains("Brain", case=False, na=False)
        else:
            mask = gtex["_primary_site"].astype(str).str.lower().eq(site.lower())
        gtex_samples_by_cancer[code] = gtex.loc[mask, "sample"].tolist()

    tcga_samples_by_cancer = {
        code: tcga.loc[tcga["cancer_code"] == code, "sample"].tolist()
        for code in sorted(set(TCGA_DISEASE_TO_CODE.values()))
    }

    return tcga_samples_by_cancer, gtex_samples_by_cancer

def load_expression_header(expr_path: Path):
    # Read only header to get sample list
    df0 = pd.read_csv(expr_path, sep="\t", nrows=0, compression="infer")
    cols = df0.columns.tolist()
    # first column is 'sample' (gene id)
    sample_cols = cols[1:]
    return cols[0], sample_cols

def compute_gene_cancer_means(expr_path: Path, tcga_samples_by_cancer, gtex_samples_by_cancer, chunksize: int = 500):
    """
    Streaming computation over the huge expression matrix.
    Expression matrix format:
      first column: gene_id (ENSG....)
      remaining columns: sample IDs
    Values are (likely) log2(TPM + offset); we'll treat them as already log-scale.
    """
    gene_col, all_sample_cols = load_expression_header(expr_path)

    # Union all needed samples
    needed_samples = set()
    for d in (tcga_samples_by_cancer, gtex_samples_by_cancer):
        for s_list in d.values():
            needed_samples.update(s_list)

    # Keep only samples present in expression matrix
    needed_samples = [s for s in all_sample_cols if s in needed_samples]

    usecols = [gene_col] + needed_samples

    # Precompute index arrays for fast mean
    sample_index = {s: i for i, s in enumerate(needed_samples)}
    tcga_idx = {code: [sample_index[s] for s in tcga_samples_by_cancer[code] if s in sample_index] for code in tcga_samples_by_cancer}
    gtex_idx = {code: [sample_index[s] for s in gtex_samples_by_cancer[code] if s in sample_index] for code in gtex_samples_by_cancer}

    rows = []  # accumulate results as list of dicts

    for chunk in pd.read_csv(expr_path, sep="\t", usecols=usecols, compression="infer", chunksize=chunksize):
        gene_ids = chunk[gene_col].astype(str).tolist()
        X = chunk[needed_samples].to_numpy(dtype=np.float32)

        for code in sorted(tcga_idx.keys()):
            t_idx = tcga_idx[code]
            n_idx = gtex_idx.get(code, [])

            # Skip if either group is empty
            if len(t_idx) == 0 or len(n_idx) == 0:
                continue

            tumor_mean = X[:, t_idx].mean(axis=1)
            normal_mean = X[:, n_idx].mean(axis=1)
            log2fc = tumor_mean - normal_mean

            # store results for this cancer for all genes in this chunk
            for g, tm, nm, fc in zip(gene_ids, tumor_mean, normal_mean, log2fc):
                rows.append({
                    "ensembl_gene": g,
                    C.cancer: code,
                    "tumor_log2tpm_mean": float(tm),
                    "normal_log2tpm_mean": float(nm),
                    C.rna_log2fc: float(fc),
                })

    out = pd.DataFrame(rows)
    return out

def build_rna_table():
    pheno = load_pheno()
    tcga_samples_by_cancer, gtex_samples_by_cancer = select_samples(pheno)
    rna = compute_gene_cancer_means(EXPR_FILE, tcga_samples_by_cancer, gtex_samples_by_cancer)
    return rna