# src/schema.py
from dataclasses import dataclass
from typing import List

@dataclass(frozen=True)
class Columns:
    # identifiers
    gene: str = "gene"
    cancer: str = "cancer"

    # RNA features (gene x cancer)
    tumor_tpm_mean: str = "tumor_tpm_mean"
    normal_tpm_mean: str = "normal_tpm_mean"
    rna_log2fc: str = "rna_log2fc"
    rna_pval: str = "rna_pval"  # optional

    # Protein target + metadata (gene x cancer)
    hpa_protein_level: str = "hpa_protein_level"  # High/Medium/Low/Not detected
    protein_detected: str = "protein_detected"    # 1/0 binary target
    hpa_n_samples: str = "hpa_n_samples"          # optional if available
    hpa_antibody_reliability: str = "hpa_antibody_reliability"  # optional

    # Gene-level features (gene only; replicated across cancers)
    gene_length: str = "gene_length"
    protein_length: str = "protein_length"
    gc_content: str = "gc_content"
    n_isoforms: str = "n_isoforms"
    tdark: str = "tdark"  # 1/0

def required_columns_for_v0() -> List[str]:
    """
    Minimal columns needed for the very first baseline model:
    predict protein_detected from rna_log2fc.
    """
    C = Columns()
    return [C.gene, C.cancer, C.rna_log2fc, C.protein_detected]