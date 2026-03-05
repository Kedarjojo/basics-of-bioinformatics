import pandas as pd

PATH = "data/raw/hpa/subcellular_location.tsv"

df = pd.read_csv(PATH, sep="\t")

print("Columns:", list(df.columns))
print("\nHead:")
print(df.head(5).to_string(index=False))

print("\nNon-null counts:")
print(df.notna().sum().sort_values(ascending=False).head(20).to_string())