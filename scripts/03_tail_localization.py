#!/usr/bin/env python3
"""
03_tail_localization.py
-----------------------
Analyzes GroEL protein sequence alignment and localizes coding sequence length differences
between mesophiles and thermophiles:
1. Calculates per-column gap frequencies and excess deletion rates (thermophile minus mesophile).
2. Maps alignment columns to E. coli GroEL numbering (548 amino acids).
3. Quantifies excess deletion mass concentrated in the intrinsically disordered C-terminal tail (>= residue 524 / column >= 540) vs folded core.
4. Computes protein length comparisons for GroEL full-length, GroEL tail, and DnaK control.

Reference:
Figure 3 and Results section of "Chaperone mRNA secondary structure does not strongly track
thermal adaptation in bacteria, but groL coding sequence shortens in its disordered C-terminal tail"
"""

import sys
import re
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy import stats

# Ensure UTF-8 output encoding across platforms
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
SEQ_DIR = BASE_DIR / "sequences"
TABLES_DIR = BASE_DIR / "tables"
TABLES_DIR.mkdir(exist_ok=True)

TAIL_COL_START = 539  # 0-indexed column 539 = alignment column 540

def clean(x):
    x = re.sub(r"[^A-Za-z0-9]+", "_", str(x).strip())
    return re.sub(r"_+", "_", x).strip("_")

def cliffs_delta(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    gt = np.sum(a[:, None] > b[None, :])
    lt = np.sum(a[:, None] < b[None, :])
    return float((gt - lt) / (len(a) * len(b)))

def analyze_tail_excess():
    meta = pd.read_csv(DATA_DIR / "groL_folding_descriptors.csv")
    type_map = dict(zip(meta.species.map(clean), meta.thermal_type))
    
    meso, thermo = [], []
    for rec in SeqIO.parse(SEQ_DIR / "GroEL_protein_aligned.fasta", "fasta"):
        tip = clean(rec.id.split()[0])
        t = type_map.get(tip)
        if t is None:
            for k, v in type_map.items():
                if tip in k or k in tip:
                    t = v
                    break
        if t == "Mesophile":
            meso.append(str(rec.seq))
        elif t == "Thermophile":
            thermo.append(str(rec.seq))
            
    aln_len = len(meso[0])
    mg = np.zeros(aln_len)
    tg = np.zeros(aln_len)
    for s in meso:
        for i, aa in enumerate(s):
            if aa == "-":
                mg[i] += 1
    for s in thermo:
        for i, aa in enumerate(s):
            if aa == "-":
                tg[i] += 1
                
    excess = tg / len(thermo) - mg / len(meso)
    pos = excess.clip(min=0)
    tot = pos.sum()
    tail = pos[TAIL_COL_START:].sum()
    pct = 100 * tail / tot
    high_excess_cols = int((excess > 0.15).sum())

    print("=" * 80)
    print("GroEL PROTEIN ALIGNMENT & C-TERMINAL TAIL LOCALIZATION:")
    print("=" * 80)
    print(f"Alignment columns: {aln_len}")
    print(f"Strains analyzed: {len(meso)} mesophiles, {len(thermo)} thermophiles")
    print(f"Total positive excess-gap mass across alignment: {tot:.3f} column-equivalents")
    print(f"Excess-gap mass in C-terminal tail (columns >= 540): {tail:.3f} column-equivalents ({pct:.1f}%)")
    print(f"Alignment columns with >15% excess deletion in thermophiles: {high_excess_cols}")
    print("-" * 80)

    # Protein length statistics table
    s4_path = DATA_DIR / "protein_length_and_tail.csv"
    if not s4_path.exists():
        s4_path = TABLES_DIR / "Table_S4_protein_length_and_tail.csv"
        
    s4 = pd.read_csv(s4_path)
    summary_rows = []
    
    print("\nPROTEIN LENGTH COMPARISONS (Mesophiles vs Thermophiles):")
    print(f"{'Protein Region':<28} {'Mesophile (mean±SD)':<22} {'Thermophile (mean±SD)':<22} {'P-value':<12} {'Cliff delta':<12}")
    print("-" * 96)
    
    for col, label in (
        ("groel_protein_length_aa", "GroEL full length"),
        ("groel_cterm_tail_aa", "GroEL C-terminal tail"),
        ("dnak_protein_length_aa", "DnaK full length (control)"),
    ):
        a = s4[s4.thermal_type == "Mesophile"][col].dropna()
        b = s4[s4.thermal_type == "Thermophile"][col].dropna()
        _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        delta = cliffs_delta(a, b)
        summary_rows.append({
            "region": label,
            "meso_mean": a.mean(),
            "meso_sd": a.std(),
            "thermo_mean": b.mean(),
            "thermo_sd": b.std(),
            "p_value": p,
            "cliffs_delta": delta
        })
        print(f"{label:<28} {a.mean():.1f} ± {a.std():.1f} aa{'':<7} {b.mean():.1f} ± {b.std():.1f} aa{'':<7} {p:<12.4g} {delta:<+12.2f}")
    print("-" * 96)

    res_df = pd.DataFrame(summary_rows)
    out_table = TABLES_DIR / "Table2_protein_lengths_tail.csv"
    res_df.to_csv(out_table, index=False)
    print(f"Saved summary table to {out_table}\n")

    return tot, tail, pct

if __name__ == "__main__":
    analyze_tail_excess()
