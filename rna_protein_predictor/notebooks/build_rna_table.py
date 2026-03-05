from src.load_rna import build_rna_table
from src.config import DATA_PROCESSED

rna = build_rna_table()

out_path = DATA_PROCESSED / "rna_gene_cancer_table.csv"
rna.to_csv(out_path, index=False)

print("Saved:", out_path)
print("Rows:", len(rna))
print("Cancers:", rna["cancer"].nunique())
print(rna["cancer"].value_counts())
print(rna.head())