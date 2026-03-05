# src/build_gene_features.py
import pandas as pd
from pathlib import Path

from .config import DATA_RAW, DATA_PROCESSED

GTF_PATH = DATA_RAW / "gene_features" / "gencode.v23.annotation.gtf"  # <-- no .gz
OUT_PATH = DATA_PROCESSED / "gene_basic_features.csv"

def parse_attributes(attr_string: str) -> dict:
    attrs = {}
    for field in attr_string.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        key, value = field.split(" ", 1)
        attrs[key] = value.replace('"', '').strip()
    return attrs

def main():
    gene_info = {}

    with open(GTF_PATH, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = parts

            if feature != "gene":
                continue

            attrs = parse_attributes(attributes)

            gene_id = attrs.get("gene_id", "").split(".")[0]  # ENSG... base
            gene_name = attrs.get("gene_name", "")
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", ""))

            if not gene_id or not gene_name:
                continue

            start_i = int(start)
            end_i = int(end)
            length = end_i - start_i + 1

            gene_info[gene_id] = {
                "gene": gene_name.upper(),
                "ensembl_base": gene_id,
                "gene_length": length,
                "gene_type": gene_type,
            }

    df = pd.DataFrame(gene_info.values())
    df.to_csv(OUT_PATH, index=False)

    print("Saved:", OUT_PATH)
    print("Genes:", len(df))
    print(df.head(5).to_string(index=False))

if __name__ == "__main__":
    main()