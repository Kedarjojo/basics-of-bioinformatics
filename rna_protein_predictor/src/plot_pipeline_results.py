import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path("data/processed")

def plot_localization_enrichment():

    df = pd.read_csv(DATA / "discordance_localization_enrichment.csv")

    overall = df[df["cancer"] == "ALL"].copy()

    overall = overall.sort_values(
        "enrichment_ratio",
        ascending=False
    ).head(15)

    plt.figure(figsize=(6,5))

    plt.barh(
        overall["feature"].str.replace("loc_",""),
        overall["enrichment_ratio"]
    )

    plt.gca().invert_yaxis()

    plt.xlabel("Discordant / Concordant enrichment")
    plt.title("Localization enrichment in RNA–protein discordance")

    plt.tight_layout()
    plt.show()


def plot_unexpected_missing_localization():

    df = pd.read_csv(
        DATA / "unexpected_missing_localization_counts.csv"
    )

    df = df.sort_values(ascending=False).head(10)

    plt.figure(figsize=(6,5))

    plt.barh(
        df.index.str.replace("loc_",""),
        df.values
    )

    plt.gca().invert_yaxis()

    plt.xlabel("Gene count")
    plt.title("Localization of unexpectedly missing proteins")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_localization_enrichment()
    plot_unexpected_missing_localization()