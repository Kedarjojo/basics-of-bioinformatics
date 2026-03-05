import pandas as pd

pheno = pd.read_csv(
    "data/raw/rnaseq/TcgaTargetGTEX_phenotype.txt",
    sep="\t",
    encoding="latin-1"   # or "ISO-8859-1"
)

print(pheno.columns)
print("\nStudy counts:")
print(pheno["_study"].value_counts())
print("\nSample type counts:")
print(pheno["_sample_type"].value_counts())
print("\nPrimary disease examples:")
print(pheno["primary disease or tissue"].unique()[:20])