#!/usr/bin/env python3
"""
plot_graphical_abstract.py
--------------------------
Generates a publication-grade, vector-based Graphical Abstract (16:9 widescreen)
summarizing the complete study in three clean panels:
1. Panel A: mRNA Secondary Structure Thermodynamics (dnaK & groL, 37°C vs 55°C, Invariant).
2. Panel B: Gene Length vs Optimal Growth Temperature (groL shortening vs dnaK control).
3. Panel C: GroEL Structural Localization (Folded core 1-523 vs disordered C-terminal tail 524-548, 73.5% deletion).

Outputs:
figures/Graphical_Abstract_Vector.png (300 DPI)
figures/Graphical_Abstract_Vector.pdf (Vector)
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Arrow
from matplotlib.lines import Line2D
from scipy import stats

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
FIG_DIR = BASE_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

def generate_vector_graphical_abstract():
    # Load data
    grol = pd.read_csv(DATA_DIR / "groL_folding_descriptors.csv")
    dnak = pd.read_csv(DATA_DIR / "dnaK_folding_descriptors.csv")
    gap = pd.read_csv(DATA_DIR / "groL_column_gap_profile.csv")

    fig = plt.figure(figsize=(16, 9), facecolor="white")
    
    # Title Banner
    fig.text(0.5, 0.94, "Chaperone mRNA Secondary Structure & GroEL C-Terminal Tail Shortening",
             fontsize=17, fontweight="bold", ha="center", color="#1A252C", fontfamily="sans-serif")
    fig.text(0.5, 0.91, "Thermal Adaptation in Bacterial Chaperones Acts on Disordered Protein Extensions Rather Than Transcript Folding",
             fontsize=11.5, fontstyle="italic", ha="center", color="#576574", fontfamily="sans-serif")

    # Create 3 sub-axes
    # [left, bottom, width, height]
    ax1 = fig.add_axes([0.05, 0.12, 0.27, 0.72])
    ax2 = fig.add_axes([0.37, 0.12, 0.27, 0.72])
    ax3 = fig.add_axes([0.69, 0.12, 0.27, 0.72])

    # Connective workflow arrows
    fig.text(0.335, 0.48, "➔", fontsize=28, ha="center", va="center", color="#BDC3C7")
    fig.text(0.655, 0.48, "➔", fontsize=28, ha="center", va="center", color="#BDC3C7")

    # =========================================================================
    # PANEL 1: mRNA Secondary Structure Invariance
    # =========================================================================
    ax1.set_title("1. mRNA Secondary Structure", fontsize=13, fontweight="bold", color="#2C3E50", pad=12)
    
    # Sub-panels for dnaK and groL MFE
    genes = ["dnaK", "groL"]
    x_pos = np.array([0, 1.2])
    width = 0.35

    # Plot normalized MFE per nucleotide at 37°C and 55°C
    meso_grol = grol[grol.thermal_type == "Mesophile"]["norm_mfe_37"].mean()
    thermo_grol = grol[grol.thermal_type == "Thermophile"]["norm_mfe_55"].mean()
    meso_dnak = dnak[dnak.thermal_type == "Mesophile"]["norm_mfe_37"].mean()
    thermo_dnak = dnak[dnak.thermal_type == "Thermophile"]["norm_mfe_55"].mean()

    m_vals = [meso_dnak, meso_grol]
    t_vals = [thermo_dnak, thermo_grol]

    bars1 = ax1.bar(x_pos - width/2, m_vals, width, label="Mesophile (37°C)", color="#3498DB", alpha=0.9, edgecolor="#2980B9")
    bars2 = ax1.bar(x_pos + width/2, t_vals, width, label="Thermophile (55°C)", color="#E74C3C", alpha=0.9, edgecolor="#C0392B")

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(["dnaK mRNA", "groL mRNA"], fontsize=11, fontweight="bold")
    ax1.set_ylabel("Normalized MFE (kcal / mol / nt)", fontsize=10.5)
    ax1.set_ylim(-0.40, 0.0)
    ax1.axhline(0, color="black", lw=0.8)
    ax1.legend(loc="lower left", fontsize=9, framealpha=0.9)

    # Invariant stamp / banner
    box_props = dict(boxstyle="round,pad=0.5", facecolor="#F8F9FA", edgecolor="#BDC3C7", lw=1.2)
    ax1.text(0.5, 0.78, "34 Statistical Tests:\nNo significant group difference\nin raw secondary structure\n(all MFE p ≥ 0.38; Cliff's δ ≈ 0)",
             transform=ax1.transAxes, ha="center", va="center", fontsize=9.5, bbox=box_props, color="#2C3E50")

    # Invariant badge
    ax1.text(0.5, 0.52, "≈ INVARIANT FOLDING ≈", transform=ax1.transAxes,
             ha="center", va="center", fontsize=10, fontweight="bold", color="#7F8C8D")

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # =========================================================================
    # PANEL 2: Continuous Length Shortening with Topt
    # =========================================================================
    ax2.set_title("2. Coding Sequence Shortening", fontsize=13, fontweight="bold", color="#2C3E50", pad=12)

    x_t = grol["topt_c"].values
    y_g = grol["cds_length_nt"].values
    y_k = dnak["cds_length_nt"].values

    # Scatter points
    ax2.scatter(x_t, y_g, color="#C0392B", s=45, alpha=0.85, label="groL CDS (shorter in thermophiles)", zorder=3)
    
    # Regression line for groL
    xs = np.linspace(x_t.min() - 2, x_t.max() + 2, 100)
    # PGLS fit line
    ax2.plot(xs, 1646.6052 - 0.3564 * xs, color="#C0392B", lw=2.4,
             label="PGLS: slope = -0.356 nt/°C\n(p = 6.0 × 10⁻⁵, λ = 0.997)", zorder=4)

    # Inset / secondary axis for dnaK control
    ax2.set_xlabel("Optimal Growth Temperature, Topt (°C)", fontsize=10.5)
    ax2.set_ylabel("groL CDS Length (nt)", fontsize=10.5, color="#C0392B")
    ax2.tick_params(axis='y', labelcolor="#C0392B")

    # Callout box for contrast
    stat_box = dict(boxstyle="round,pad=0.5", facecolor="#FDEDEC", edgecolor="#E6B0AA", lw=1.2)
    ax2.text(0.5, 0.20, "groL: 1636.3 nt ➔ 1626.5 nt\nCliff's δ = 0.571 (p = 5.1 × 10⁻⁴)\nContinuous across all 4 phyla\n\ndnaK Control: Invariant (p = 0.71)",
             transform=ax2.transAxes, ha="center", va="center", fontsize=9.5, bbox=stat_box, color="#922B21")

    ax2.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
    ax2.spines["top"].set_visible(False)
    ax2.grid(True, linestyle=":", alpha=0.4)

    # =========================================================================
    # PANEL 3: GroEL Tail Localization Mechanism
    # =========================================================================
    ax3.set_title("3. GroEL C-Terminal Tail Truncation", fontsize=13, fontweight="bold", color="#2C3E50", pad=12)

    # Gap frequency profile
    diff_col = "diff_thermo_minus_meso"
    col_idx = "column"
    x_cols = gap[col_idx].values
    y_diff = gap[diff_col].values

    # Color bars: slate for core (<540), gold for tail (>=540)
    colors = ["#F39C12" if c >= 540 else "#7F8C8D" for c in x_cols]
    ax3.bar(x_cols, y_diff, color=colors, width=1.1, zorder=3)
    ax3.axhline(0, color="black", lw=0.8, zorder=2)
    ax3.axvline(540, color="#D4AC0D", linestyle="--", lw=1.5, zorder=4)

    ax3.set_xlim(1, len(gap) + 5)
    ax3.set_ylim(-0.15, 0.40)
    ax3.set_xlabel("GroEL Alignment Column (MAFFT)", fontsize=10.5)
    ax3.set_ylabel("Excess Deletion in Thermophiles\n(Thermophile − Mesophile Gap %)", fontsize=10.5)

    # Highlight shaded regions
    ax3.axvspan(1, 539, color="#ECF0F1", alpha=0.5, zorder=0)
    ax3.axvspan(540, len(gap) + 5, color="#FEF9E7", alpha=0.7, zorder=0)

    # Annotations
    ax3.text(270, 0.34, "Folded Core (Residues 1–523)\nIntact & Conserved", ha="center", fontsize=9.5, fontweight="bold", color="#7F8C8D")
    ax3.text(555, 0.28, "C-Tail (524–548)\n73.5% Deletion", ha="center", fontsize=9.5, fontweight="bold", color="#B7950B")

    tail_box = dict(boxstyle="round,pad=0.5", facecolor="#FEF9E7", edgecolor="#F9E79F", lw=1.2)
    ax3.text(0.5, 0.18, "Tail Length:\nMeso: 21.8 ± 3.3 aa  vs  Thermo: 18.9 ± 2.3 aa\n(p = 1.7 × 10⁻³, Cliff's δ = +0.52)\nShortening of client-filtering GGM tail",
             transform=ax3.transAxes, ha="center", va="center", fontsize=9.5, bbox=tail_box, color="#7D6608")

    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.grid(axis="y", linestyle=":", alpha=0.4)

    # Bottom summary footer
    fig.text(0.5, 0.03, "Summary: Bacterial chaperone mRNA structure does not track thermal niche; instead, GroEL selectively shortens its disordered C-terminal peptide tail.",
             fontsize=10.5, fontweight="bold", ha="center", color="#2C3E50", fontfamily="sans-serif")

    out_png = FIG_DIR / "Graphical_Abstract_Vector.png"
    out_pdf = FIG_DIR / "Graphical_Abstract_Vector.pdf"
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()
    print(f"Generated Vector Graphical Abstract: {out_png} and {out_pdf}")

if __name__ == "__main__":
    generate_vector_graphical_abstract()
