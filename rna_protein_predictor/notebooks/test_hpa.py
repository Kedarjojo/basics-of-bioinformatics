from src.load_hpa import load_hpa_clean

df = load_hpa_clean()

print(df.head())
print("\nCounts:")
print(df["protein_detected"].value_counts())
print("\nCancers:")
print(df["cancer"].unique())