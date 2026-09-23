#!/usr/bin/env python3
"""
plot_figure3.py
---------------
Generates Figure 3 of the manuscript:
Per-column difference in gap frequency (thermophile minus mesophile) across the GroEL protein alignment.
- Panel A: Full protein mapped to E. coli residue numbering (1-548). Shaded regions denote
  the folded core (residues 1-523, grey) vs the intrinsically disordered C-terminal tail (residues 524-548, gold).
- Panel B: High-resolution zoom of alignment columns >= 500 displaying the concentration of
  excess deletion signal in the C-terminal tail (~73.5% of total signal).

Outputs:
figures/Fig3_tail_localisation_v3_final.png
figures/Fig3_tail_localisation_v3_final.pdf
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
FIG_DIR = BASE_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

def generate_figure3():
    gp_path = DATA_DIR / "groL_column_gap_profile.csv"
    gp = pd.read_csv(gp_path)
    
    # Column names handling
    res_col = "ecoli_resnum" if "ecoli_resnum" in gp.columns else "ecoli_residue"
    diff_col = "diff_thermo_minus_meso" if "diff_thermo_minus_meso" in gp.columns else "gap_freq_diff"
    col_idx = "column" if "column" in gp.columns else "alignment_column"

    gp["x"] = gp[res_col].interpolate(limit_direction="both")
    L = len(gp)
    C524 = 539  # 0-indexed column 539 corresponds to column 540 / residue 524
    pos = gp[diff_col].clip(lower=0)
    tot, tail = pos.sum(), pos[C524:].sum()

    fig = plt.figure(figsize=(13, 5.5), constrained_layout=True)
    gs = GridSpec(1, 2, width_ratios=[2, 1], figure=fig)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])

    # ---------- Panel A: Full Alignment ----------
    axA.axvspan(1, 523, color="#eef1f4", zorder=0, label="Folded core (1-523)")
    axA.axvspan(524, 548, color="#ffe08a", alpha=0.55, zorder=0, label="C-terminal tail (524-548)")
    axA.bar(gp["x"], gp[diff_col], width=1.6, color="#b5342a", linewidth=0, zorder=3)
    axA.axhline(0, color="black", linewidth=0.9, zorder=2)

    axA.set_xlim(1, 552)
    axA.set_ylim(-0.20, 0.40)
    axA.set_xlabel("E. coli GroEL residue position", fontsize=11)
    axA.set_ylabel("Gap-frequency difference\n(thermophile − mesophile)", fontsize=11)
    axA.set_title("A   Folded core (grey) vs. C-terminal tail (gold)", loc="left", fontsize=11.5, fontweight="bold")
    for s in ("top", "right"):
        axA.spines[s].set_visible(False)
    axA.grid(axis="y", linestyle=":", alpha=0.4)

    # ---------- Panel B: Zoom of Terminal Region ----------
    sub = gp.iloc[500:].copy()
    axB.axvspan(500, C524, color="#eef1f4", zorder=0)
    axB.axvspan(C524, L + 1, color="#ffe08a", alpha=0.55, zorder=0)
    colors = ["#b5342a" if i >= C524 else "#9aa7b4" for i in sub.index]
    x_sub = sub[col_idx] if col_idx in sub.columns else sub.index + 1
    axB.bar(x_sub, sub[diff_col], width=0.92, color=colors, linewidth=0, zorder=3)
    axB.axhline(0, color="black", linewidth=0.9, zorder=2)
    axB.axvline(C524 + 0.5, color="#8a6100", linestyle="--", linewidth=1.2, zorder=4)

    axB.set_xlim(500, L + 1)
    axB.set_ylim(-0.08, 0.40)
    axB.set_xlabel("Alignment column", fontsize=11)
    axB.set_title(f"B   Tail: {tail:.2f}/{tot:.2f} = {100*tail/tot:.0f}% of signal", loc="left", fontsize=11.5, fontweight="bold")
    for s in ("top", "right"):
        axB.spines[s].set_visible(False)
    axB.grid(axis="y", linestyle=":", alpha=0.4)

    out_png = FIG_DIR / "Fig3_tail_localisation_v3_final.png"
    out_pdf = FIG_DIR / "Fig3_tail_localisation_v3_final.pdf"
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()
    print(f"Generated Figure 3: {out_png} and {out_pdf} (tail excess: {tail:.2f}/{tot:.2f} = {100*tail/tot:.1f}%)")

if __name__ == "__main__":
    generate_figure3()
