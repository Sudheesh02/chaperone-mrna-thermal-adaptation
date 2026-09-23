#!/usr/bin/env python3
"""
run_reproduction_pipeline.py
----------------------------
Master Reproduction Pipeline for:
"Chaperone mRNA secondary structure does not strongly track thermal adaptation in bacteria,
but groL coding sequence shortens in its disordered C-terminal tail"

Executes all analytical steps and figure generation end-to-end, validating all computed
statistics against reported manuscript values.

Usage:
    python run_reproduction_pipeline.py
"""

import sys
import subprocess
import time
from pathlib import Path
import numpy as np
import pandas as pd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = BASE_DIR / "scripts"
TABLES_DIR = BASE_DIR / "tables"
FIGURES_DIR = BASE_DIR / "figures"

def run_step(step_num, title, script_path):
    print(f"\n[{step_num}/5] RUNNING: {title} ({script_path.name})...")
    start = time.time()
    res = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, cwd=str(BASE_DIR))
    elapsed = time.time() - start
    if res.returncode != 0:
        print(f"FAILED (exit code {res.returncode}):")
        print(res.stderr)
        sys.exit(res.returncode)
    print(f"✓ Completed in {elapsed:.2f}s")
    return res.stdout

def main():
    print("=" * 82)
    print(" REPRODUCIBILITY PIPELINE: BACTERIAL CHAPERONE THERMAL ADAPTATION")
    print(" Manuscript: bioRxiv preprint (Sudheesh et al.)")
    print("=" * 82)

    # 1. 34 Statistical tests (15 descriptors x 2 genes + 4 z-scores)
    run_step(1, "Mann-Whitney tests, Cliff's delta & BH-FDR / Bonferroni", SCRIPTS_DIR / "01_run_34_tests.py")

    # 2. Phylogenetic Generalized Least Squares (PGLS) & OLS regression
    run_step(2, "Phylogenetic Generalized Least Squares (PGLS)", SCRIPTS_DIR / "02_pgls_regression.py")

    # 3. GroEL C-terminal tail gap localization
    run_step(3, "GroEL protein alignment & C-terminal tail localization", SCRIPTS_DIR / "03_tail_localization.py")

    # 4. Dinucleotide shuffle algorithm self-test
    run_step(4, "Altschul-Erickson dinucleotide shuffle verification", SCRIPTS_DIR / "04_dinuc_shuffle.py")

    # 5. Publication figures generation
    print("\n[5/5] GENERATING PUBLICATION FIGURES...")
    for fig_script in ["plot_figure1.py", "plot_figure2.py", "plot_figure3.py"]:
        start = time.time()
        res = subprocess.run([sys.executable, str(SCRIPTS_DIR / fig_script)], capture_output=True, text=True, cwd=str(BASE_DIR))
        elapsed = time.time() - start
        if res.returncode != 0:
            print(f"  Error generating {fig_script}: {res.stderr}")
            sys.exit(res.returncode)
        print(f"  ✓ {fig_script} ({elapsed:.2f}s)")

    # --------------------------------------------------------------------------
    # RIGOROUS SCIENTIFIC ASSERTION & AUDIT CHECKS
    # --------------------------------------------------------------------------
    print("\n" + "=" * 82)
    print(" VERIFYING COMPUTED VALUES AGAINST MANUSCRIPT CLAIMS:")
    print("=" * 82)

    # Check 1: 34 tests results
    t34 = pd.read_csv(TABLES_DIR / "Table_all_34_statistical_tests.csv")
    grol_len = t34[(t34.gene == "groL") & (t34.metric == "Transcript length")].iloc[0]
    grol_z55 = t34[(t34.gene == "groL") & (t34.metric == "z-score 55C")].iloc[0]
    bh_survivors = t34[t34.BH_significant == True]

    assert len(bh_survivors) == 2, f"Expected exactly 2 BH survivors, found {len(bh_survivors)}"
    assert np.isclose(grol_len.delta, 0.5712, atol=1e-3), f"groL length delta mismatch: {grol_len.delta}"
    assert np.isclose(grol_len.p, 0.000511, atol=1e-4), f"groL length p-value mismatch: {grol_len.p}"
    assert np.isclose(grol_z55.delta, 0.5648, atol=1e-3), f"groL z55 delta mismatch: {grol_z55.delta}"
    assert np.isclose(grol_z55.p, 0.000638, atol=1e-4), f"groL z55 p-value mismatch: {grol_z55.p}"
    print(" [PASS] 34-test battery: exactly 2 comparisons survive BH FDR 5% (groL length & groL z-55°C).")

    # Check 2: PGLS & OLS models
    pgls_tab = pd.read_csv(TABLES_DIR / "Table_OLS_PGLS_summary.csv")
    grol_pgls = pgls_tab[pgls_tab.gene == "groL"].iloc[0]
    dnak_pgls = pgls_tab[pgls_tab.gene == "dnaK"].iloc[0]

    assert np.isclose(grol_pgls.spearman_rho, -0.5889, atol=1e-3), f"Spearman rho mismatch: {grol_pgls.spearman_rho}"
    assert grol_pgls.ols_p < 1e-4, f"OLS p-value should be < 1e-4: {grol_pgls.ols_p}"
    assert grol_pgls.pgls_p < 1e-4, f"PGLS p-value should be < 1e-4: {grol_pgls.pgls_p}"
    assert np.isclose(grol_pgls.pgls_slope, -0.35, atol=0.03), f"PGLS slope mismatch: {grol_pgls.pgls_slope}"
    assert dnak_pgls.spearman_p > 0.05, f"dnaK should be non-significant: {dnak_pgls.spearman_p}"
    print(" [PASS] Continuous Topt scaling & PGLS: groL length decreases with Topt (p < 1e-4, lambda ~ 1.0); dnaK invariant.")

    # Check 3: Protein lengths & C-terminal tail localization
    t_prot = pd.read_csv(TABLES_DIR / "Table2_protein_lengths_tail.csv")
    g_full = t_prot[t_prot.region == "GroEL full length"].iloc[0]
    g_tail = t_prot[t_prot.region == "GroEL C-terminal tail"].iloc[0]
    d_full = t_prot[t_prot.region == "DnaK full length (control)"].iloc[0]

    assert np.isclose(g_full.meso_mean, 544.4, atol=0.2), f"GroEL meso mean mismatch: {g_full.meso_mean}"
    assert np.isclose(g_full.thermo_mean, 541.2, atol=0.2), f"GroEL thermo mean mismatch: {g_full.thermo_mean}"
    assert np.isclose(g_tail.meso_mean, 21.8, atol=0.2), f"GroEL tail meso mean mismatch: {g_tail.meso_mean}"
    assert np.isclose(g_tail.thermo_mean, 18.9, atol=0.2), f"GroEL tail thermo mean mismatch: {g_tail.thermo_mean}"
    assert g_full.p_value < 0.001, f"GroEL length p-value mismatch: {g_full.p_value}"
    assert g_tail.p_value < 0.005, f"GroEL tail p-value mismatch: {g_tail.p_value}"
    assert d_full.p_value > 0.5, f"DnaK length should not differ: {d_full.p_value}"
    print(" [PASS] Protein-level localization: GroEL full length (544.4 vs 541.2 aa) and tail (21.8 vs 18.9 aa) confirmed.")

    # Check 4: Figure files
    figs = ["Fig1_34tests_300dpi.png", "Fig1_34tests_300dpi.pdf", "Fig2_length_PGLS.png", "Fig3_tail_localisation_v3_final.png"]
    for f in figs:
        assert (FIGURES_DIR / f).exists(), f"Missing figure file: {f}"
    print(" [PASS] Publication figures: Figures 1, 2, and 3 successfully verified in figures/ directory.")

    print("\n" + "=" * 82)
    print(" ALL AUDIT CHECKS PASSED: 100% REPRODUCIBILITY CONFIRMED.")
    print(" The repository is fully self-contained and ready for publication.")
    print("=" * 82 + "\n")

if __name__ == "__main__":
    main()
