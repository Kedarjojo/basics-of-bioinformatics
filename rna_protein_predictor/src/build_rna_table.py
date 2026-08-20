"""Command-line entry point for the streamed Xena RNA summary."""

from .config import DATA_PROCESSED
from .load_rna import build_rna_table


def main() -> None:
    output = DATA_PROCESSED / "rna_gene_cancer_table.csv"
    table = build_rna_table()
    table.to_csv(output, index=False)
    print(f"Saved {output} ({len(table):,} rows)")


if __name__ == "__main__":
    main()

