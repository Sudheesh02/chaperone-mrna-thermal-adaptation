#!/usr/bin/env python3
"""
plot_figure1.py
---------------
Generates Figure 1 of the manuscript:
Effect sizes (Cliff's delta: mesophile - thermophile) and two-sided Mann-Whitney p-values
for all 34 pre-specified comparisons between mesophiles and thermophiles, displayed
separately for dnaK (left) and groL (right).
Red indicates tests surviving Benjamini-Hochberg FDR correction at alpha = 0.05.
Grey indicates non-surviving tests.

Outputs:
figures/Fig1_34tests_300dpi.png
figures/Fig1_34tests_300dpi.pdf
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
TABLES_DIR = BASE_DIR / "tables"
FIG_DIR = BASE_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

def generate_figure1():
    table_path = TABLES_DIR / "Table_all_34_statistical_tests.csv"
    if not table_path.exists():
        # Fall back to Table_S2 if not yet run
        table_path = TABLES_DIR / "Table_S2_statistics_34tests.csv"
    
    res = pd.read_csv(table_path)
    # Ensure BH_significant is boolean
    res["BH_significant"] = res["BH_significant"].astype(bool)

    fig, axes = plt.subplots(1, 2, figsize=(12, 9), sharex=True)

    for ax, gene in zip(axes, ["dnaK", "groL"]):
        d = res[res.gene == gene].sort_values("delta").reset_index(drop=True)
        y = np.arange(len(d))

        for i, row in d.iterrows():
            color = "#c0392b" if row.BH_significant else "#7f8c8d"
            ax.plot(row.delta, i, "o", color=color, markersize=7.5, zorder=3)
            # Annotate p-value on right side
            p_val = row.p if "p" in row else row["p_value"]
            ax.text(0.72, i, f"p={p_val:.3f}", fontsize=7.5, va="center", color="#444444")

        ax.axvline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(d.metric, fontsize=8.5)
        ax.set_xlabel("Cliff’s δ (mesophile − thermophile)", fontsize=10.5)
        ax.set_xlim(-0.75, 0.90)
        ax.set_title(gene, fontsize=14, fontweight="bold", style="italic", pad=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="x", linestyle=":", alpha=0.4)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#c0392b", markersize=9, label="Survives BH-FDR 5%"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#7f8c8d", markersize=9, label="Does not survive"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=2, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, -0.01))

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.08)
    
    out_png = FIG_DIR / "Fig1_34tests_300dpi.png"
    out_pdf = FIG_DIR / "Fig1_34tests_300dpi.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Generated Figure 1: {out_png} and {out_pdf}")

if __name__ == "__main__":
    generate_figure1()
