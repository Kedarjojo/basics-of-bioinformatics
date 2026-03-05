import pandas as pd

from .config import DATA_PROCESSED

FILE = DATA_PROCESSED / "unexpected_missing_proteins_top.csv"

def main():

    df = pd.read_csv(FILE)

    loc_cols = [c for c in df.columns if c.startswith("loc_")]

    print("\nLocalization frequencies:")
    loc_counts = df[loc_cols].sum().sort_values(ascending=False)
    print(loc_counts.head(15))

    print("\nMean expression:", df["tumor_log2tpm_mean"].mean())
    print("Median expression:", df["tumor_log2tpm_mean"].median())

    print("\nTop cancers:")
    print(df["cancer"].value_counts())

    # export localization summary
    loc_counts.to_csv(
        DATA_PROCESSED / "unexpected_missing_localization_counts.csv"
    )

    print("\nSaved localization summary")

if __name__ == "__main__":
    main()