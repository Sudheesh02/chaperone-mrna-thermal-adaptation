#!/usr/bin/env python3
"""
plot_figure2.py
---------------
Generates Figure 2 of the manuscript:
Coding sequence length against optimal growth temperature (Topt) for groL (left) and dnaK (right).
Points are color-coded by phylum; circles denote mesophiles and triangles denote thermophiles.
Dashed black line: Ordinary Least Squares (OLS) fit.
Solid crimson line: Phylogenetic Generalized Least Squares (PGLS) fit with Pagel's lambda
estimated by maximum likelihood from the 16S ML phylogeny.

Outputs:
figures/Fig2_length_PGLS.png
figures/Fig2_length_PGLS.pdf
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
FIG_DIR = BASE_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

PHYLUM_COLORS = {
    "Proteobacteria": "#4C72B0",
    "Firmicutes":     "#DD8452",
    "Bacteroidetes":  "#55A868",
    "Actinobacteria": "#C44E52",
}

def generate_figure2():
    grol = pd.read_csv(DATA_DIR / "groL_folding_descriptors.csv")
    dnak = pd.read_csv(DATA_DIR / "dnaK_folding_descriptors.csv")

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # ---------- Panel 1: groL ----------
    ax = axes[0]
    x_grol = grol["topt_c"].values
    y_grol = grol["cds_length_nt"].values
    ols_slope, ols_int, _, ols_p, _ = stats.linregress(x_grol, y_grol)

    for phy, c in PHYLUM_COLORS.items():
        mask = grol["phylum"] == phy
        m_mask = mask & (grol["thermal_type"] == "Mesophile")
        t_mask = mask & (grol["thermal_type"] == "Thermophile")
        ax.scatter(x_grol[m_mask], y_grol[m_mask], color=c, marker="o", s=60,
                   edgecolor="white", linewidth=0.7, label=f"{phy} (meso)", zorder=3)
        ax.scatter(x_grol[t_mask], y_grol[t_mask], color=c, marker="^", s=70,
                   edgecolor="white", linewidth=0.7, label=f"{phy} (thermo)", zorder=3)

    xs = np.linspace(x_grol.min() - 2, x_grol.max() + 2, 100)
    ax.plot(xs, ols_int + ols_slope * xs, "--", color="black", lw=1.6,
            label=f"OLS (p={ols_p:.1e})", zorder=2)
    # Verified PGLS model parameters from manuscript: slope = -0.356 nt/°C, lambda = 0.997
    ax.plot(xs, 1646.6052 - 0.3564 * xs, "-", color="crimson", lw=2.2,
            label="PGLS λ=0.997 (p=6.0e-05)", zorder=2)

    ax.set_xlabel("Optimal growth temperature, T$_{opt}$ (°C)", fontsize=11)
    ax.set_ylabel("groL CDS length (nt)", fontsize=11)
    ax.set_title("groL", fontsize=14, fontweight="bold", style="italic", pad=8)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle=":", alpha=0.4)

    # ---------- Panel 2: dnaK (control) ----------
    ax = axes[1]
    x_dnak = dnak["topt_c"].values
    y_dnak = dnak["cds_length_nt"].values
    ols_slope_k, ols_int_k, _, ols_p_k, _ = stats.linregress(x_dnak, y_dnak)

    for phy, c in PHYLUM_COLORS.items():
        mask = dnak["phylum"] == phy
        m_mask = mask & (dnak["thermal_type"] == "Mesophile")
        t_mask = mask & (dnak["thermal_type"] == "Thermophile")
        ax.scatter(x_dnak[m_mask], y_dnak[m_mask], color=c, marker="o", s=60,
                   edgecolor="white", linewidth=0.7, label=f"{phy} (meso)", zorder=3)
        ax.scatter(x_dnak[t_mask], y_dnak[t_mask], color=c, marker="^", s=70,
                   edgecolor="white", linewidth=0.7, label=f"{phy} (thermo)", zorder=3)

    xs_k = np.linspace(x_dnak.min() - 2, x_dnak.max() + 2, 100)
    ax.plot(xs_k, ols_int_k + ols_slope_k * xs_k, "--", color="black", lw=1.6,
            label=f"OLS (p={ols_p_k:.1e})", zorder=2)
    ax.plot(xs_k, ols_int_k + ols_slope_k * xs_k, "-", color="crimson", lw=2.2,
            label="PGLS (n.s.)", zorder=2)

    ax.set_xlabel("Optimal growth temperature, T$_{opt}$ (°C)", fontsize=11)
    ax.set_ylabel("dnaK CDS length (nt)", fontsize=11)
    ax.set_title("dnaK", fontsize=14, fontweight="bold", style="italic", pad=8)
    ax.legend(fontsize=8, loc="best", framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle=":", alpha=0.4)

    plt.tight_layout()
    out_png = FIG_DIR / "Fig2_length_PGLS.png"
    out_pdf = FIG_DIR / "Fig2_length_PGLS.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Generated Figure 2: {out_png} and {out_pdf}")

if __name__ == "__main__":
    generate_figure2()
