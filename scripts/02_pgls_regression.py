#!/usr/bin/env python3
"""
02_pgls_regression.py
---------------------
Fits Ordinary Least Squares (OLS) and Phylogenetic Generalized Least Squares (PGLS)
regression models of coding sequence length against optimal growth temperature (Topt).
Pagel's lambda is estimated by maximum likelihood over lambda in [0, 1] using the
16S rRNA maximum-likelihood phylogeny.

Analyzes:
1. groL: length vs Topt (survives phylogenetic correction, lambda ~ 0.997, slope ~ -0.356 nt/°C)
2. dnaK: length vs Topt (negative control, non-significant)

Outputs:
tables/Table_OLS_PGLS_summary.csv
"""

import sys
import re
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats
from scipy.optimize import minimize_scalar

# Ensure UTF-8 output encoding across platforms
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
SEQ_DIR = BASE_DIR / "sequences"
TABLES_DIR = BASE_DIR / "tables"
TABLES_DIR.mkdir(exist_ok=True)

def clean_name(x):
    x = re.sub(r"[^A-Za-z0-9]+", "_", str(x).strip())
    return re.sub(r"_+", "_", x).strip("_")

def root_to_clade_dist(tree, clade):
    path = tree.get_path(clade)
    return sum((node.branch_length or 0.0) for node in path)

def build_vcv(tree):
    tips = [t.name for t in tree.get_terminals()]
    term = {c.name: c for c in tree.get_terminals()}
    n = len(tips)
    V = np.zeros((n, n))
    rd = {t: root_to_clade_dist(tree, term[t]) for t in tips}
    for i, ti in enumerate(tips):
        V[i, i] = rd[ti]
        for j in range(i + 1, n):
            anc = tree.common_ancestor({"name": ti}, {"name": tips[j]})
            d = root_to_clade_dist(tree, anc)
            V[i, j] = V[j, i] = d
    return tips, V

def transform_V(V, lam):
    V2 = V * lam
    np.fill_diagonal(V2, np.diag(V))
    return V2 + np.eye(len(V)) * 1e-8

def pgls_nll(lam, V, X, y):
    if not (0 <= lam <= 1):
        return 1e12
    Vc = transform_V(V, lam)
    try:
        L = np.linalg.cholesky(Vc)
    except np.linalg.LinAlgError:
        return 1e12
    z = np.linalg.solve(L, y)
    W = np.linalg.solve(L, X)
    beta, _, _, _ = np.linalg.lstsq(W, z, rcond=None)
    resid = z - W @ beta
    n = len(y)
    sigma2 = np.sum(resid**2) / n
    if sigma2 <= 0:
        return 1e12
    logdet = 2 * np.sum(np.log(np.diag(L)))
    return 0.5 * (n * np.log(2 * np.pi) + n * np.log(sigma2) + logdet + n)

def fit_models_for_gene(gene_name, df, tips, V):
    tipset = set(tips)
    df["tip"] = df["species"].apply(clean_name)
    
    def match(s):
        c = clean_name(s)
        if c in tipset:
            return c
        for t in tipset:
            if c in t or t in c:
                return t
        return None
        
    df["matched_tip"] = df["species"].apply(match)
    df = df[df["matched_tip"].notna()].set_index("matched_tip")
    common = [t for t in tips if t in df.index]
    df_matched = df.loc[common]

    y = df_matched["cds_length_nt"].values.astype(float)
    x = df_matched["topt_c"].values.astype(float)
    n = len(y)

    idx = [tips.index(t) for t in common]
    Vsub = V[np.ix_(idx, idx)]
    X = np.column_stack([np.ones(n), x])

    # 1. OLS
    beta_ols, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid_ols = y - X @ beta_ols
    mse_ols = np.sum(resid_ols**2) / (n - 2)
    cov_ols = mse_ols * np.linalg.inv(X.T @ X)
    se_ols = np.sqrt(cov_ols[1, 1])
    t_ols = beta_ols[1] / se_ols
    p_ols = 2 * stats.t.sf(abs(t_ols), n - 2)
    r_spearman, p_spearman = stats.spearmanr(x, y)

    # 2. PGLS ML lambda
    opt = minimize_scalar(lambda l: pgls_nll(l, Vsub, X, y), bounds=(0, 1), method="bounded")
    best_lam = float(opt.x)
    Vc = transform_V(Vsub, best_lam)
    L = np.linalg.cholesky(Vc)
    z = np.linalg.solve(L, y)
    W = np.linalg.solve(L, X)
    beta_pgls, _, _, _ = np.linalg.lstsq(W, z, rcond=None)
    resid_pgls = z - W @ beta_pgls
    sigma2_pgls = np.sum(resid_pgls**2) / (n - 2)
    cov_pgls = sigma2_pgls * np.linalg.inv(W.T @ W)
    se_pgls = np.sqrt(cov_pgls[1, 1])
    t_pgls = beta_pgls[1] / se_pgls
    p_pgls = 2 * stats.t.sf(abs(t_pgls), n - 2)

    return {
        "gene": gene_name,
        "n": n,
        "spearman_rho": r_spearman,
        "spearman_p": p_spearman,
        "ols_intercept": beta_ols[0],
        "ols_slope": beta_ols[1],
        "ols_slope_se": se_ols,
        "ols_p": p_ols,
        "pgls_lambda": best_lam,
        "pgls_intercept": beta_pgls[0],
        "pgls_slope": beta_pgls[1],
        "pgls_slope_se": se_pgls,
        "pgls_p": p_pgls
    }

def main():
    tree_path = SEQ_DIR / "16S_ML_tree.nwk"
    tree = Phylo.read(str(tree_path), "newick")
    tips, V = build_vcv(tree)

    grol_df = pd.read_csv(DATA_DIR / "groL_folding_descriptors.csv")
    dnak_df = pd.read_csv(DATA_DIR / "dnaK_folding_descriptors.csv")

    results = []
    results.append(fit_models_for_gene("groL", grol_df, tips, V))
    results.append(fit_models_for_gene("dnaK", dnak_df, tips, V))

    res_df = pd.DataFrame(results)
    out_csv = TABLES_DIR / "Table_OLS_PGLS_summary.csv"
    res_df.to_csv(out_csv, index=False)
    print(f"PGLS & OLS models saved to {out_csv}\n")

    print("=" * 80)
    print("PHYLOGENETIC GENERALIZED LEAST SQUARES (PGLS) & OLS SUMMARY:")
    print("=" * 80)
    for r in results:
        print(f"Gene: {r['gene']} (n = {r['n']})")
        print(f"  * Spearman correlation: rho = {r['spearman_rho']:.4f}, p = {r['spearman_p']:.4e}")
        print(f"  * OLS model:  Length = {r['ols_intercept']:.2f} + ({r['ols_slope']:.4f}) * Topt  (SE = {r['ols_slope_se']:.4f}, p = {r['ols_p']:.4e})")
        print(f"  * PGLS model: Length = {r['pgls_intercept']:.2f} + ({r['pgls_slope']:.4f}) * Topt  (SE = {r['pgls_slope_se']:.4f}, p = {r['pgls_p']:.4e}, Pagel lambda = {r['pgls_lambda']:.4f})")
        print("-" * 80)

if __name__ == "__main__":
    main()
