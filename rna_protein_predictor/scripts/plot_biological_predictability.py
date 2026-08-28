#!/usr/bin/env python3
"""Create publication figures for CPTAC biological predictability results."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"

BLUE = "#2471A3"
LIGHT_BLUE = "#85C1E9"
ORANGE = "#D97706"
RED = "#B03A2E"
GRAY = "#6B7280"
LIGHT_GRAY = "#D1D5DB"
BLACK = "#1F2937"

MODEL_LABELS = {
    "intercept_only": "Intercept only",
    "biological": "Biological only",
    "technical": "Technical only",
    "combined": "Technical + biological",
}

COMPARISON_LABELS = {
    "full_vs_technical": "All biology added",
    "reduced_vs_technical": "Reduced-overlap biology added",
    "localization_group_vs_technical": "Localization added",
    "unique_localization_full_vs_no_localization": "Unique localization contribution",
    "length_group_vs_technical": "Protein length added",
    "unique_length_full_vs_no_length": "Unique length contribution",
    "topology_group_vs_technical": "Topology added",
    "unique_topology_full_vs_no_topology": "Unique topology contribution",
    "full_vs_reduced_overlap": "Broad overlapping terms added",
}

FEATURE_LABELS = {
    "log1p_protein_length": "Protein length (log₁₊)",
    "log1p_transmembrane_count": "Transmembrane count (log₁₊)",
    "has_signal_peptide": "Signal peptide",
    "uniprot_secreted": "Secreted",
    "uniprot_extracellular": "Extracellular",
    "uniprot_plasma_membrane": "Plasma membrane",
    "uniprot_membrane_any": "Any membrane",
    "uniprot_cytosol": "Cytosol",
    "uniprot_nucleus": "Nucleus",
    "uniprot_mitochondrion": "Mitochondrion",
    "uniprot_endoplasmic_reticulum": "Endoplasmic reticulum",
    "uniprot_golgi": "Golgi",
    "uniprot_lysosome": "Lysosome",
    "uniprot_peroxisome": "Peroxisome",
    "uniprot_cell_junction": "Cell junction",
}


def configure_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "axes.edgecolor": BLACK,
        "axes.linewidth": 0.8,
        "xtick.color": BLACK,
        "ytick.color": BLACK,
        "text.color": BLACK,
        "axes.labelcolor": BLACK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def require(path: Path, columns: set[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing input table: {path}")
    frame = pd.read_csv(path)
    missing = sorted(columns.difference(frame.columns))
    if missing:
        raise ValueError(f"{path.name} is missing columns: {missing}")
    return frame


def clean_axis(ax: plt.Axes, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(axis="x", color=LIGHT_GRAY, linewidth=0.6, alpha=0.65)
        ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")


def save(fig: plt.Figure, stem: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight",
                facecolor="white")


def model_performance_panel(ax: plt.Axes, intervals: pd.DataFrame) -> None:
    data = intervals.loc[
        (intervals["metric"] == "r2") &
        intervals["model"].isin(["biological", "technical", "combined"])
    ].copy()
    order = ["biological", "technical", "combined"]
    data["order"] = data["model"].map({name: i for i, name in enumerate(order)})
    data = data.sort_values("order")
    y = np.arange(len(data))
    colors = [LIGHT_BLUE, GRAY, BLUE]
    errors = np.vstack([
        data["estimate"] - data["ci_low"],
        data["ci_high"] - data["estimate"],
    ])
    ax.errorbar(data["estimate"], y, xerr=errors, fmt="none", ecolor=BLACK,
                elinewidth=1.2, capsize=3, zorder=2)
    ax.scatter(data["estimate"], y, s=65, c=colors, edgecolor="white",
               linewidth=0.8, zorder=3)
    for yi, value in zip(y, data["estimate"]):
        ax.text(value + 0.012, yi, f"{value:.3f}", va="center", fontsize=8.5)
    ax.set_yticks(y, [MODEL_LABELS[x] for x in data["model"]])
    ax.set_xlim(0, max(0.44, data["ci_high"].max() + 0.045))
    ax.set_xlabel("Held-out gene R²")
    ax.set_title("Biology improves held-out prediction", loc="left")
    ax.invert_yaxis()
    clean_axis(ax)


def ablation_panel(ax: plt.Axes, paired: pd.DataFrame) -> None:
    wanted = list(COMPARISON_LABELS)
    data = paired.loc[
        (paired["metric"] == "r2") & paired["comparison"].isin(wanted)
    ].copy()
    data["order"] = data["comparison"].map({name: i for i, name in enumerate(wanted)})
    data = data.sort_values("order")
    y = np.arange(len(data))
    errors = np.vstack([
        data["estimate"] - data["ci_low"],
        data["ci_high"] - data["estimate"],
    ])
    colors = [
        BLUE if name in {"full_vs_technical", "reduced_vs_technical"}
        else GRAY if name == "full_vs_reduced_overlap"
        else ORANGE
        for name in data["comparison"]
    ]
    ax.axvline(0, color=BLACK, linewidth=0.9, linestyle="--")
    ax.errorbar(data["estimate"], y, xerr=errors, fmt="none", ecolor=BLACK,
                elinewidth=1.0, capsize=2.5, zorder=2)
    ax.scatter(data["estimate"], y, s=48, c=colors, edgecolor="white",
               linewidth=0.7, zorder=3)
    ax.set_yticks(y, [COMPARISON_LABELS[x] for x in data["comparison"]])
    ax.set_xlabel("Favorable change in held-out R²")
    ax.set_title("Localization provides the largest biological increment", loc="left")
    ax.invert_yaxis()
    clean_axis(ax)


def coefficient_panel(ax: plt.Axes, stability: pd.DataFrame) -> None:
    biological = set(FEATURE_LABELS)
    data = stability.loc[
        (stability["model"] == "combined") & stability["feature"].isin(biological)
    ].copy().sort_values("median_standardized_coefficient")
    y = np.arange(len(data))
    errors = np.vstack([
        data["median_standardized_coefficient"] - data["ci_low"],
        data["ci_high"] - data["median_standardized_coefficient"],
    ])
    colors = np.where(data["median_standardized_coefficient"] >= 0, BLUE, RED)
    ax.axvline(0, color=BLACK, linewidth=0.9, linestyle="--")
    ax.errorbar(data["median_standardized_coefficient"], y, xerr=errors,
                fmt="none", ecolor=BLACK, elinewidth=0.9, capsize=2, zorder=2)
    ax.scatter(data["median_standardized_coefficient"], y, c=colors, s=42,
               edgecolor="white", linewidth=0.6, zorder=3)
    ax.set_yticks(y, [FEATURE_LABELS[x] for x in data["feature"]])
    ax.set_xlabel("Standardized ridge coefficient")
    ax.set_title("Adjusted biological coefficient stability", loc="left")
    clean_axis(ax)


def descriptive_panel(ax: plt.Axes, groups: pd.DataFrame) -> None:
    data = groups.copy().sort_values("present_minus_absent_median")
    y = np.arange(len(data))
    errors = np.vstack([
        data["present_minus_absent_median"] - data["bootstrap_ci_low"],
        data["bootstrap_ci_high"] - data["present_minus_absent_median"],
    ])
    colors = np.where(data["present_minus_absent_median"] >= 0, BLUE, RED)
    ax.axvline(0, color=BLACK, linewidth=0.9, linestyle="--")
    ax.errorbar(data["present_minus_absent_median"], y, xerr=errors,
                fmt="none", ecolor=BLACK, elinewidth=0.9, capsize=2, zorder=2)
    ax.scatter(data["present_minus_absent_median"], y, c=colors, s=42,
               edgecolor="white", linewidth=0.6, zorder=3)
    labels = [f"{FEATURE_LABELS[x]} (n={n:,})"
              for x, n in zip(data["feature"], data["present_n"])]
    ax.set_yticks(y, labels)
    ax.set_xlabel("Median predictability difference: annotated − other genes")
    ax.set_title("Unadjusted localization contrasts", loc="left")
    clean_axis(ax)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=FIGURES)
    args = parser.parse_args()
    configure_style()

    model_intervals = require(
        TABLES / "gene_predictability_biology_model_intervals.csv",
        {"model", "metric", "estimate", "ci_low", "ci_high"},
    )
    paired = require(
        TABLES / "gene_predictability_biology_sensitivity_paired_intervals.csv",
        {"comparison", "metric", "estimate", "ci_low", "ci_high"},
    )
    stability = require(
        TABLES / "gene_predictability_biology_coefficient_stability.csv",
        {"model", "feature", "median_standardized_coefficient", "ci_low", "ci_high"},
    )
    groups = require(
        TABLES / "gene_predictability_biology_group_summary.csv",
        {"feature", "present_n", "present_minus_absent_median",
         "bootstrap_ci_low", "bootstrap_ci_high"},
    )

    individual = [
        ("biology_model_performance", model_performance_panel, model_intervals, (6.4, 3.5)),
        ("biology_grouped_ablation", ablation_panel, paired, (8.2, 5.0)),
        ("biology_coefficient_stability", coefficient_panel, stability, (7.8, 6.0)),
        ("biology_localization_descriptive", descriptive_panel, groups, (8.2, 5.8)),
    ]
    for stem, function, data, size in individual:
        fig, ax = plt.subplots(figsize=size)
        function(ax, data)
        fig.tight_layout()
        save(fig, stem, args.output_dir)
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 11.5))
    model_performance_panel(axes[0, 0], model_intervals)
    ablation_panel(axes[0, 1], paired)
    coefficient_panel(axes[1, 0], stability)
    descriptive_panel(axes[1, 1], groups)
    for ax, label in zip(axes.flat, "ABCD"):
        panel_label(ax, label)
    fig.suptitle(
        "Technical and biological determinants of RNA–protein predictability",
        fontsize=15, fontweight="bold", x=0.5, y=1.005,
    )
    fig.text(
        0.5, -0.012,
        "n = 8,672 genes. Model intervals use 2,000 paired gene bootstraps; "
        "coefficient ranges show the 2.5th–97.5th percentiles across 100 CV fits. "
        "Panel D is unadjusted.",
        ha="center", va="bottom", fontsize=8, color=GRAY,
    )
    fig.tight_layout(h_pad=2.2, w_pad=3.0)
    save(fig, "biology_summary_composite", args.output_dir)
    plt.close(fig)

    print(f"Figures written to {args.output_dir.resolve()}")
    for stem, *_ in individual:
        print(f"  {stem}.png / {stem}.pdf")
    print("  biology_summary_composite.png / biology_summary_composite.pdf")


if __name__ == "__main__":
    main()
