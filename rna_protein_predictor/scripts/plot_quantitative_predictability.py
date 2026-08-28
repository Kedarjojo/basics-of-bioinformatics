#!/usr/bin/env python3
"""Create publication figures for the CPTAC quantitative benchmark.

The script intentionally reads the analysis tables rather than recomputing any
models.  It therefore preserves the frozen patient folds, bootstrap draws, and
the gene-level atlas used by the primary analyses.
"""

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
PURPLE = "#6C3483"
GRAY = "#6B7280"
LIGHT_GRAY = "#D1D5DB"
BLACK = "#1F2937"
CANCERS = ["BRCA", "COAD", "GBM", "PDAC"]
MODEL_COLORS = {"ols": BLUE, "ridge": PURPLE, "huber": ORANGE}
MODEL_LABELS = {"ols": "OLS (primary)", "ridge": "Ridge", "huber": "Huber"}


def configure_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.labelsize": 9.5, "axes.edgecolor": BLACK, "axes.linewidth": 0.8,
        "xtick.color": BLACK, "ytick.color": BLACK, "text.color": BLACK,
        "axes.labelcolor": BLACK, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def require(name: str, columns: set[str]) -> pd.DataFrame:
    path = TABLES / name
    if not path.exists():
        raise FileNotFoundError(
            f"Missing input table: {path}. Run the corresponding analysis stage first."
        )
    data = pd.read_csv(path)
    missing = sorted(columns.difference(data.columns))
    if missing:
        raise ValueError(f"{name} is missing columns: {missing}")
    return data


def clean_axis(ax: plt.Axes, grid: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=LIGHT_GRAY, linewidth=0.6, alpha=0.7)
        ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.13, 1.06, label, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")


def save(fig: plt.Figure, stem: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")


def bootstrap_panel(ax: plt.Axes, intervals: pd.DataFrame, estimand: str,
                    title: str, xlabel: str, percent: bool = True) -> None:
    data = intervals.loc[(intervals["bootstrap"] == "patient") &
                         (intervals["model"] == "ols") &
                         (intervals["estimand"] == estimand)].copy()
    data["cancer"] = pd.Categorical(data["cancer"], CANCERS, ordered=True)
    data = data.sort_values("cancer")
    y = np.arange(len(data))
    errors = np.vstack([data["estimate"] - data["ci_low"], data["ci_high"] - data["estimate"]])
    ax.errorbar(data["estimate"], y, xerr=errors, fmt="none", ecolor=BLACK,
                elinewidth=1.1, capsize=3, zorder=2)
    ax.scatter(data["estimate"], y, s=62, color=BLUE, edgecolor="white", linewidth=0.8, zorder=3)
    for yi, value in zip(y, data["estimate"]):
        label = f"{value:.1%}" if percent else f"{value:.3f}"
        ax.annotate(label, (value, yi), xytext=(7, 0), textcoords="offset points", va="center", fontsize=8)
    ax.set_yticks(y, data["cancer"])
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left")
    clean_axis(ax, "x")


def sensitivity_panel(ax: plt.Axes, summary: pd.DataFrame) -> None:
    data = summary.copy()
    data["cancer"] = pd.Categorical(data["cancer"], CANCERS, ordered=True)
    for model in ("ols", "ridge", "huber"):
        d = data.loc[data["model"] == model].sort_values("cancer")
        ax.plot(d["cancer"], d["relative_mse_improvement"], marker="o", markersize=5,
                linewidth=1.6, color=MODEL_COLORS[model], label=MODEL_LABELS[model])
    ax.set_ylim(0.12, 0.36)
    ax.set_ylabel("Relative MSE improvement")
    ax.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    ax.set_title("Estimator choice does not alter the conclusion", loc="left")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    clean_axis(ax)


def gene_bootstrap_panel(ax: plt.Axes, intervals: pd.DataFrame, estimand: str,
                         title: str, xlabel: str) -> None:
    data = intervals.loc[(intervals["bootstrap"] == "gene") &
                         (intervals["model"] == "ols") &
                         (intervals["estimand"] == estimand)].copy()
    data["cancer"] = pd.Categorical(data["cancer"], CANCERS, ordered=True)
    data = data.sort_values("cancer")
    y = np.arange(len(data))
    errors = np.vstack([data["estimate"] - data["ci_low"], data["ci_high"] - data["estimate"]])
    ax.errorbar(data["estimate"], y, xerr=errors, fmt="none", ecolor=BLACK,
                elinewidth=1.1, capsize=3, zorder=2)
    ax.scatter(data["estimate"], y, s=58, color=ORANGE, edgecolor="white", linewidth=0.8, zorder=3)
    ax.axvline(0, color=BLACK, linestyle="--", linewidth=0.8)
    ax.set_yticks(y, data["cancer"])
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.xaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    ax.set_title(title, loc="left")
    clean_axis(ax, "x")


def cancer_performance_figure(summary: pd.DataFrame, intervals: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.0))
    bootstrap_panel(axes[0, 0], intervals, "relative_mse_improvement",
                    "RNA improves held-out protein prediction", "Patient-bootstrap relative MSE improvement")
    sensitivity_panel(axes[0, 1], summary)
    gene_bootstrap_panel(axes[1, 0], intervals, "median_gene_improvement",
                         "Typical gene-level benefit is positive", "Gene-bootstrap median MSE improvement")
    gene_bootstrap_panel(axes[1, 1], intervals, "fraction_genes_improved",
                         "Most genes improve over their training-fold mean", "Gene-bootstrap fraction of genes improved")
    for ax, label in zip(axes.flat, "ABCD"):
        panel_label(ax, label)
    fig.suptitle("RNA predicts quantitative protein abundance in all four CPTAC cancers",
                 fontsize=14.5, fontweight="bold", y=1.01)
    fig.text(0.5, -0.015, "Five-fold held-out-patient evaluation; OLS is primary. "
             "Bars show 95% percentile intervals from 2,000 patient (A) or gene (C–D) bootstraps.",
             ha="center", fontsize=8, color=GRAY)
    fig.tight_layout(h_pad=2.5, w_pad=3.0)
    save(fig, "quantitative_cancer_performance", out)
    plt.close(fig)


def stability_heatmap(ax: plt.Axes, pairs: pd.DataFrame, title: str) -> None:
    matrix = pd.DataFrame(np.eye(4), index=CANCERS, columns=CANCERS)
    for row in pairs.itertuples(index=False):
        matrix.loc[row.cancer_a, row.cancer_b] = row.gene_profile_spearman
        matrix.loc[row.cancer_b, row.cancer_a] = row.gene_profile_spearman
    image = ax.imshow(matrix, cmap="Blues", vmin=0, vmax=1)
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{matrix.iloc[i, j]:.2f}", ha="center", va="center",
                    color="white" if matrix.iloc[i, j] > 0.58 else BLACK, fontsize=9)
    ax.set_xticks(range(4), CANCERS)
    ax.set_yticks(range(4), CANCERS)
    ax.set_title(title, loc="left")
    plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Spearman correlation")


def loco_distribution_panel(ax: plt.Axes, loco: pd.DataFrame) -> None:
    ordered = [loco.loc[loco["held_out_cancer"] == cancer, "test_spearman"].dropna().to_numpy() for cancer in CANCERS]
    violin = ax.violinplot(ordered, showmedians=True, showextrema=False)
    for body in violin["bodies"]:
        body.set_facecolor(LIGHT_BLUE); body.set_edgecolor(BLUE); body.set_alpha(0.85)
    violin["cmedians"].set_color(BLACK)
    ax.axhline(0, color=BLACK, linestyle="--", linewidth=0.8)
    ax.set_xticks(range(1, 5), CANCERS)
    ax.set_ylabel("Gene-level test Spearman correlation")
    ax.set_title("Coupling direction transfers to unseen cancers", loc="left")
    clean_axis(ax)


def loco_positive_panel(ax: plt.Axes, loco: pd.DataFrame) -> None:
    data = loco.groupby("held_out_cancer", as_index=False)["test_spearman"].agg(
        fraction_positive=lambda x: (x > 0).mean(), median="median", genes="count")
    data["held_out_cancer"] = pd.Categorical(data["held_out_cancer"], CANCERS, ordered=True)
    data = data.sort_values("held_out_cancer")
    bars = ax.bar(data["held_out_cancer"], data["fraction_positive"], color=BLUE, width=0.65)
    for bar, value in zip(bars, data["fraction_positive"]):
        ax.text(bar.get_x() + bar.get_width()/2, value + .012, f"{value:.1%}", ha="center", fontsize=8.5)
    ax.set_ylim(0.85, 1.0)
    ax.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    ax.set_ylabel("Genes with positive transfer correlation")
    ax.set_title("Positive direction is retained for >92% of genes", loc="left")
    clean_axis(ax)


def transfer_figure(pairs: pd.DataFrame, loco: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1, 1.25, 1]})
    stability_heatmap(axes[0], pairs, "Observed coupling profiles are concordant")
    loco_distribution_panel(axes[1], loco)
    loco_positive_panel(axes[2], loco)
    for ax, label in zip(axes, "ABC"):
        panel_label(ax, label)
    fig.suptitle("Gene-level RNA–protein coupling generalizes across cancers", fontsize=14.5, fontweight="bold", y=1.03)
    fig.text(0.5, -0.04, "Panel A: pairwise cross-cancer concordance of observed gene-level RNA–protein Spearman correlations. "
             "Panels B–C: leave-one-cancer-out standardized slopes, assessed by within-gene patient ranking in the held-out cancer.",
             ha="center", fontsize=8, color=GRAY)
    fig.tight_layout(w_pad=2.7)
    save(fig, "quantitative_cross_cancer_transfer", out)
    plt.close(fig)


def atlas_distribution_panel(ax: plt.Axes, atlas: pd.DataFrame) -> None:
    for count, color, label in [(3, LIGHT_BLUE, "Observed in 3 cancers"), (4, BLUE, "Observed in 4 cancers")]:
        values = atlas.loc[atlas["cancers_observed"] == count, "primary_median_ridge_improvement"]
        ax.hist(values, bins=np.linspace(-0.30, 0.90, 49), density=True, alpha=0.72, color=color, label=label)
    ax.axvline(0, color=BLACK, linestyle="--", linewidth=0.8)
    ax.set_xlabel("Median cross-cancer Ridge MSE improvement")
    ax.xaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    ax.set_ylabel("Gene density")
    ax.set_title("Most genes show positive cross-cancer predictability", loc="left")
    ax.legend(frameon=False, fontsize=8)
    clean_axis(ax)


def per_cancer_distribution_panel(ax: plt.Axes, by_cancer: pd.DataFrame) -> None:
    # Relative MSE can have rare, arbitrarily large negative values when a
    # gene's baseline error is very small.  Show robust distribution summaries
    # rather than allowing a few ratios to collapse the central 95% of genes.
    rows = []
    for cancer in CANCERS:
        values = by_cancer.loc[
            by_cancer["cancer"] == cancer, "ridge_mse_improvement"
        ].dropna()
        rows.append({
            "cancer": cancer,
            "median": values.median(),
            "q25": values.quantile(0.25), "q75": values.quantile(0.75),
            "q05": values.quantile(0.05), "q95": values.quantile(0.95),
        })
    data = pd.DataFrame(rows)
    x = np.arange(len(data))
    ax.vlines(x, data["q05"], data["q95"], color=LIGHT_BLUE, linewidth=2.0, zorder=1)
    ax.vlines(x, data["q25"], data["q75"], color=BLUE, linewidth=7.5, zorder=2)
    ax.scatter(x, data["median"], color=BLACK, s=34, zorder=3)
    ax.axhline(0, color=BLACK, linestyle="--", linewidth=0.8)
    ax.set_xticks(x, CANCERS)
    ax.set_ylabel("Gene-level Ridge MSE improvement")
    ax.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    ax.set_title("Cohort variation is not driven by extreme ratios", loc="left")
    clean_axis(ax)


def atlas_rank_stability_panel(ax: plt.Axes, by_cancer: pd.DataFrame) -> None:
    wide = by_cancer.pivot(index="gene", columns="cancer", values="ridge_mse_improvement").reindex(columns=CANCERS)
    matrix = wide.corr(method="spearman", min_periods=100)
    image = ax.imshow(matrix, cmap="Purples", vmin=0, vmax=1)
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{matrix.iloc[i, j]:.2f}", ha="center", va="center",
                    color="white" if matrix.iloc[i, j] > 0.55 else BLACK, fontsize=9)
    ax.set_xticks(range(4), CANCERS)
    ax.set_yticks(range(4), CANCERS)
    ax.set_title("Predictability rankings are stable across cancers", loc="left")
    plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Spearman correlation")


def directional_consistency_panel(ax: plt.Axes, atlas: pd.DataFrame) -> None:
    rows = []
    for count in (3, 4):
        d = atlas.loc[atlas["cancers_observed"] == count]
        rows.append({"label": f"{count} cancers", "all_positive": (d["positive_ridge_cancers"] == count).mean(),
                     "mixed_or_negative": (d["positive_ridge_cancers"] < count).mean()})
    data = pd.DataFrame(rows)
    ax.bar(data["label"], data["all_positive"], color=BLUE, label="Positive in every observed cancer")
    ax.bar(data["label"], data["mixed_or_negative"], bottom=data["all_positive"], color=LIGHT_GRAY,
           label="At least one non-positive cancer")
    for idx, value in enumerate(data["all_positive"]):
        ax.text(idx, value / 2, f"{value:.1%}", color="white", ha="center", va="center", fontsize=9, fontweight="bold")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    ax.set_ylabel("Genes")
    ax.set_title("Positive predictability is common across coverage strata", loc="left")
    ax.legend(frameon=False, fontsize=7.5, loc="lower left")
    clean_axis(ax)


def atlas_figure(atlas: pd.DataFrame, by_cancer: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2))
    atlas_distribution_panel(axes[0, 0], atlas)
    per_cancer_distribution_panel(axes[0, 1], by_cancer)
    atlas_rank_stability_panel(axes[1, 0], by_cancer)
    directional_consistency_panel(axes[1, 1], atlas)
    for ax, label in zip(axes.flat, "ABCD"):
        panel_label(ax, label)
    fig.suptitle("An 8,681-gene atlas reveals reproducible RNA-to-protein predictability", fontsize=14.5, fontweight="bold", y=1.01)
    fig.text(0.5, -0.012, "Primary target: median cross-cancer Ridge MSE improvement versus the training-fold gene mean. "
             "Panel B shows medians (points), IQRs (thick bars), and 5th–95th percentiles (thin bars), limiting the influence of rare unstable ratios. "
             "Genes were evaluated in at least three cancers; no high/low predictability threshold was used for modeling.",
             ha="center", fontsize=8, color=GRAY)
    fig.tight_layout(h_pad=2.5, w_pad=3.0)
    save(fig, "quantitative_gene_predictability_atlas", out)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=FIGURES)
    args = parser.parse_args()
    configure_style()

    summary = require("cptac_quantitative_model_summary.csv", {"cancer", "model", "relative_mse_improvement"})
    intervals = require("cptac_quantitative_bootstrap_intervals.csv", {"bootstrap", "cancer", "model", "estimand", "estimate", "ci_low", "ci_high"})
    pairs = require("cptac_coupling_pairwise_stability.csv", {"cancer_a", "cancer_b", "gene_profile_spearman"})
    loco = require("cptac_quantitative_loco_transfer.csv", {"held_out_cancer", "gene", "test_spearman"})
    atlas = require("cptac_gene_predictability_atlas.csv", {"gene", "cancers_observed", "primary_median_ridge_improvement", "positive_ridge_cancers"})
    by_cancer = require("cptac_gene_predictability_by_cancer.csv", {"cancer", "gene", "ridge_mse_improvement"})

    cancer_performance_figure(summary, intervals, args.output_dir)
    transfer_figure(pairs, loco, args.output_dir)
    atlas_figure(atlas, by_cancer, args.output_dir)
    print(f"Figures written to {args.output_dir.resolve()}")
    for stem in ("quantitative_cancer_performance", "quantitative_cross_cancer_transfer", "quantitative_gene_predictability_atlas"):
        print(f"  {stem}.png / {stem}.pdf")


if __name__ == "__main__":
    main()
