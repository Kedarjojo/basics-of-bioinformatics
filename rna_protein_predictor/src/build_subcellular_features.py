# src/build_subcellular_features.py
import re
import pandas as pd
from collections import Counter
from .config import DATA_RAW, DATA_PROCESSED

IN_PATH = DATA_RAW / "hpa" / "subcellular_location.tsv"
OUT_WIDE = DATA_PROCESSED / "hpa_subcellular_features_wide.csv"
OUT_LONG = DATA_PROCESSED / "hpa_subcellular_features_long.csv"

def _pick_col(cols, preferred):
    lower = {c.lower(): c for c in cols}
    for key in preferred:
        if key.lower() in lower:
            return lower[key.lower()]
    # fallback: partial match
    for c in cols:
        for key in preferred:
            if key.lower() in c.lower():
                return c
    return None

def _normalize_token(tok: str) -> str:
    tok = tok.strip()
    tok = re.sub(r"\s+", " ", tok)
    # make consistent column-safe token
    tok = tok.replace("/", " or ")
    tok = re.sub(r"[^A-Za-z0-9 _\-]", "", tok)
    tok = tok.strip().replace(" ", "_").replace("-", "_")
    return tok

def _split_locations(text: str) -> list[str]:
    if text is None:
        return []
    s = str(text).strip()
    if s == "" or s.lower() in {"nan", "none"}:
        return []
    # common separators in HPA exports
    parts = re.split(r"[;,|]", s)
    out = []
    for p in parts:
        p = p.strip()
        if p:
            out.append(p)
    return out

def main():
    df = pd.read_csv(IN_PATH, sep="\t")

    gene_col = _pick_col(df.columns, ["Gene name", "Gene", "gene", "Symbol", "GeneName"])
    main_col = "Main location"
    add_col = "Additional location"
    extra_col = "Extracellular location"

    df["_combined_locations"] = (
            df[main_col].fillna("").astype(str) + ";" +
            df[add_col].fillna("").astype(str) + ";" +
            df[extra_col].fillna("").astype(str)
    )
    loc_col = "_combined_locations"

    if gene_col is None:
        raise ValueError(f"Could not detect gene column in {IN_PATH}. Columns: {list(df.columns)}")
    if loc_col is None:
        # some exports split main/additional; try combining
        main_col = _pick_col(df.columns, ["Main location", "Main_location"])
        add_col = _pick_col(df.columns, ["Additional location", "Additional_location"])
        if main_col and add_col:
            df["_combined_locations"] = df[main_col].astype(str).fillna("") + ";" + df[add_col].astype(str).fillna("")
            loc_col = "_combined_locations"
        else:
            raise ValueError(f"Could not detect localization column in {IN_PATH}. Columns: {list(df.columns)}")

    # Normalize gene symbols
    df["_gene"] = df[gene_col].astype(str).str.upper().str.strip()

    # Build long table: one row per (gene, location)
    rows = []
    for g, txt in zip(df["_gene"], df[loc_col]):
        locs = _split_locations(txt)
        for loc in locs:
            rows.append({"gene": g, "location_raw": loc, "location": _normalize_token(loc)})

    long = pd.DataFrame(rows).drop_duplicates()
    long.to_csv(OUT_LONG, index=False)

    # keep only locations that appear often enough to be meaningful (avoid 1-off noise)
    loc_counts = Counter(long["location"].tolist())
    common_locs = {loc for loc, n in loc_counts.items() if n >= 50}

    long_common = long[long["location"].isin(common_locs)].copy()

    wide = (
        long_common.assign(value=1)
        .pivot_table(index="gene", columns="location", values="value", aggfunc="max", fill_value=0)
        .reset_index()
    )

    # prefix columns
    wide = wide.rename(columns={c: (c if c == "gene" else f"loc_{c}") for c in wide.columns})

    wide.to_csv(OUT_WIDE, index=False)

    print("Detected gene column:", gene_col)
    print("Detected location column:", loc_col)
    print("Saved long:", OUT_LONG, "| rows:", len(long))
    print("Saved wide:", OUT_WIDE, "| genes:", wide["gene"].nunique(), "| loc features:", wide.shape[1] - 1)
    print("Top locations:")
    top = pd.Series(loc_counts).sort_values(ascending=False).head(15)
    print(top.to_string())

if __name__ == "__main__":
    main()