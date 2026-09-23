#!/usr/bin/env python3
"""
01_run_34_tests.py
------------------
Recomputes the 34 pre-specified comparisons between mesophilic and thermophilic bacteria:
- 15 thermodynamic/structural descriptors per gene (dnaK, groL) = 30 primary tests
- 4 dinucleotide-controlled structural z-scores (37°C, 55°C) = 4 tests
Applies two-sided Mann-Whitney U tests, computes Cliff's delta effect sizes,
and performs Benjamini-Hochberg FDR (alpha = 0.05) and Bonferroni multiple testing corrections.

Reference:
Table 1 and Results section of "Chaperone mRNA secondary structure does not strongly track
thermal adaptation in bacteria, but groL coding sequence shortens in its disordered C-terminal tail"
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

# Ensure UTF-8 output encoding across platforms
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
TABLES_DIR = BASE_DIR / "tables"
TABLES_DIR.mkdir(exist_ok=True)

def cliffs_delta(a, b):
    """Computes Cliff's delta dominance statistic between two independent groups."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    gt = np.sum(a[:, None] > b[None, :])
    lt = np.sum(a[:, None] < b[None, :])
    return float((gt - lt) / (len(a) * len(b)))

def load_data(filename):
    return pd.read_csv(DATA_DIR / filename)

METRICS = [
    ("MFE 37C", "mfe_37_kcal_mol"),
    ("MFE 55C", "mfe_55_kcal_mol"),
    ("norm. MFE 37C", "norm_mfe_37"),
    ("norm. MFE 55C", "norm_mfe_55"),
    ("Delta-Delta-G", "ddg_mfe_55_minus_37"),
    ("GC%", "gc_pct"),
    ("%% paired 37C", "pct_paired_mfe_37"),
    ("%% paired 55C", "pct_paired_mfe_55"),
    ("Ensemble gap 37C", "ensemble_gap_37_kcal_mol"),
    ("Ensemble gap 55C", "ensemble_gap_55_kcal_mol"),
    ("MFE-centroid gap 37C", "mfe_centroid_gap_37"),
    ("MFE-centroid gap 55C", "mfe_centroid_gap_55"),
    ("Transcript length", "cds_length_nt"),
    ("Ensemble diversity 37C (per-nt)", "ed_norm_37"),
    ("Ensemble diversity 55C (per-nt)", "ed_norm_55"),
]

def run_analysis():
    rows = []
    # 1. 30 Primary thermodynamic and structural descriptors
    for gene, fname in (("dnaK", "dnaK_folding_descriptors.csv"),
                        ("groL", "groL_folding_descriptors.csv")):
        df = load_data(fname)
        meso = df[df.thermal_type == "Mesophile"]
        thermo = df[df.thermal_type == "Thermophile"]
        for label, col in METRICS:
            a = meso[col].dropna()
            b = thermo[col].dropna()
            stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            delta = cliffs_delta(a, b)
            rows.append({
                "gene": gene,
                "metric": label,
                "p": p,
                "delta": delta,
                "n_meso": len(a),
                "n_thermo": len(b)
            })

    # 2. 4 Dinucleotide-controlled structural z-scores (N = 100 shuffles)
    for gene, fname in (("groL", "groL_zscores_N100.csv"),
                        ("dnaK", "dnaK_zscores_N100.csv")):
        df = load_data(fname)
        for temp in (37, 55):
            sub = df[df.temp == temp]
            type_col = "Type" if "Type" in sub.columns else "thermal_type"
            z_col = "z" if "z" in sub.columns else "z_score"
            a = sub[sub[type_col] == "Mesophile"][z_col].dropna()
            b = sub[sub[type_col] == "Thermophile"][z_col].dropna()
            stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            delta = cliffs_delta(a, b)
            rows.append({
                "gene": gene,
                "metric": f"z-score {temp}C",
                "p": p,
                "delta": delta,
                "n_meso": len(a),
                "n_thermo": len(b)
            })

    res = pd.DataFrame(rows)
    m = len(res)  # 34 tests

    # Benjamini-Hochberg FDR correction at alpha = 0.05
    res = res.sort_values("p").reset_index(drop=True)
    res["bh_rank"] = np.arange(1, m + 1)
    res["bh_critical"] = (res["bh_rank"] / m) * 0.05
    res["bonferroni_p"] = np.minimum(res["p"] * m, 1.0)

    # BH step-up procedure
    significant_ranks = res[res["p"] <= res["bh_critical"]]["bh_rank"]
    max_k = significant_ranks.max() if not significant_ranks.empty else 0
    res["BH_significant"] = res["bh_rank"] <= max_k
    res["Bonferroni_significant"] = res["p"] <= (0.05 / m)

    out_csv = TABLES_DIR / "Table_all_34_statistical_tests.csv"
    res.to_csv(out_csv, index=False)
    print(f"Results successfully saved to {out_csv}")

    print("\n" + "=" * 80)
    print("34 STATISTICAL COMPARISONS SUMMARY (sorted by p-value):")
    print("=" * 80)
    print(f"{'Gene':<6} {'Metric':<34} {'p-value':<12} {'Cliff delta':<12} {'BH 5%':<8} {'Bonferroni':<10}")
    print("-" * 80)
    for _, r in res.iterrows():
        bh_str = "YES" if r["BH_significant"] else "No"
        bonf_str = "YES" if r["Bonferroni_significant"] else "No"
        print(f"{r['gene']:<6} {r['metric']:<34} {r['p']:<12.6g} {r['delta']:<12.4f} {bh_str:<8} {bonf_str:<10}")
    print("-" * 80)

    survivors = res[res["BH_significant"]]
    print(f"\nTests surviving Benjamini-Hochberg FDR (alpha = 0.05): {len(survivors)}")
    for _, r in survivors.iterrows():
        print(f"  * {r['gene']} - {r['metric']}: p = {r['p']:.6g}, Cliff's delta = {r['delta']:.4f}")

    return res

if __name__ == "__main__":
    run_analysis()
