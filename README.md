# Chaperone mRNA Secondary Structure & GroEL C-Terminal Tail Shortening Under Thermal Adaptation

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Preprint: bioRxiv](https://img.shields.io/badge/bioRxiv-Preprint-red.svg)](#citation)

Companion code, dataset, and reproduction repository for the short communication manuscript:

> **"Chaperone mRNA secondary structure does not strongly track thermal adaptation in bacteria, but *groL* coding sequence shortens in its disordered C-terminal tail"**  
> *Sudheesh Singh, Vaibhav Sahu, Sushree Sradha, Vedikaa Dhiman* (bioRxiv, 2026)

---

## Scientific Overview

Ribosomal and transfer RNAs show established compositional and structural signatures of adaptation to growth temperature ($T_{opt}$), but messenger RNA generally does not. In this study, we investigated whether transcripts for the two principal bacterial chaperone systems, *dnaK* (Hsp70) and *groL* (GroEL/Hsp60), encode systematically different mRNA secondary structures in thermophiles relative to mesophiles.

Using a phylum-matched cohort of 50 bacterial strains (25 mesophiles, 25 thermophiles balanced across Actinobacteria, Bacteroidetes, Firmicutes, and Proteobacteria), we folded each coding sequence at 37°C and 55°C and evaluated 15 thermodynamic and structural descriptors per gene.

### Key Findings

1. **Raw mRNA folding does not separate thermal groups**: Out of 30 primary thermodynamic descriptors, no minimum free energy (MFE) measure (raw or normalized per nucleotide) differed detectably between mesophiles and thermophiles after multiple-testing correction (all $P \ge 0.38$, Mann-Whitney $U$).
2. **Residual signal after dinucleotide control**: After controlling for dinucleotide composition with Altschul-Erickson shuffles ($N = 100$), thermophile *groL* retained modestly more residual secondary structure at 55°C (mean $z = -1.41$ vs $-2.77$; $P = 6.4 \times 10^{-4}$, Cliff's $\delta = +0.56$), surviving Benjamini-Hochberg FDR at 5%.
3. **Dominant signal is groL coding sequence shortening**: Thermophile *groL* coding sequences are significantly shorter (mean $1626.5 \pm 7.2$ nt vs $1636.3 \pm 10.1$ nt; Cliff's $\delta = 0.571$, $P = 5.1 \times 10^{-4}$). This scales continuously with $T_{opt}$ (Spearman $\rho = -0.589$, $P = 6.8 \times 10^{-6}$) and persists after phylogenetic generalized least squares correction ($n = 50$, Pagel's $\lambda = 0.997$, slope $= -0.356\text{ nt/}^\circ\text{C}$, $P = 6.0 \times 10^{-5}$ under R `nlme::gls`). The *dnaK* control shows no comparable relationship (Spearman $\rho = -0.057, P = 0.70$).
4. **Localization to the disordered C-terminal tail**: Multiple sequence alignment of all 50 GroEL proteins localizes ~73.5% of the excess deletion signal directly to the terminal ~25 residues (columns $\ge 540$, corresponding to residues 524–548 in *E. coli* numbering), which forms the intrinsically disordered, glycine/methionine-rich C-terminal tail.

---

## Repository Structure

```text
chaperone-mrna-thermal-adaptation/
├── README.md                      # Academic documentation & study guide
├── LICENSE                        # Open-source MIT License
├── requirements.txt               # Pinned Python package dependencies
├── environment.yml                # Conda environment specification
├── run_reproduction_pipeline.py   # End-to-end reproduction & validation pipeline
│
├── data/                          # Curated machine-readable datasets
│   ├── strain_metadata.csv        # 50 bacterial strains (taxonomy, Topt, Tmin, Tmax)
│   ├── groL_folding_descriptors.csv # 15 thermodynamic descriptors at 37°C and 55°C
│   ├── dnaK_folding_descriptors.csv # 15 thermodynamic descriptors at 37°C and 55°C
│   ├── groL_zscores_N100.csv      # groL shuffle z-scores (N = 100)
│   ├── groL_zscores_N50.csv       # groL shuffle z-scores (N = 50 robustness)
│   ├── groL_zscores_N20.csv       # groL shuffle z-scores (N = 20 robustness)
│   ├── dnaK_zscores_N100.csv      # dnaK shuffle z-scores (N = 100)
│   ├── dnaK_zscores_N50.csv       # dnaK shuffle z-scores (N = 50 robustness)
│   ├── dnaK_zscores_N20.csv       # dnaK shuffle z-scores (N = 20 robustness)
│   ├── groL_column_gap_profile.csv# MAFFT alignment per-column gap differences
│   └── protein_length_and_tail.csv# GroEL full length, tail length, and DnaK length
│
├── sequences/                     # Sequence files & phylogenetic tree
│   ├── groL_CDS.fasta             # 50 groL coding sequences
│   ├── dnaK_CDS.fasta             # 50 dnaK coding sequences
│   ├── GroEL_protein_unaligned.fasta
│   ├── GroEL_protein_aligned.fasta# MAFFT alignment (569 columns)
│   ├── DnaK_protein_unaligned.fasta
│   ├── 16S_rRNA_unaligned.fasta   # 16S rRNA gene sequences
│   ├── 16S_rRNA_aligned.fasta     # MAFFT aligned 16S rRNA sequences
│   └── 16S_ML_tree.nwk            # FastTree maximum-likelihood phylogeny (Newick)
│
├── scripts/                       # Modular analysis & visualization scripts
│   ├── 01_run_34_tests.py         # 34 Mann-Whitney tests, Cliff's delta & BH-FDR
│   ├── 02_pgls_regression.py      # Python PGLS & OLS models with ML Pagel's lambda
│   ├── 02_pgls_regression.R       # Companion R script (ape + nlme::gls)
│   ├── 03_tail_localization.py    # Per-column gap profiles & 73.5% tail excess
│   ├── 04_dinuc_shuffle.py        # Altschul-Erickson dinucleotide shuffle algorithm
│   ├── plot_figure1.py            # Generates Figure 1 (Cliff's delta forest plot)
│   ├── plot_figure2.py            # Generates Figure 2 (Length vs Topt + PGLS lines)
│   └── plot_figure3.py            # Generates Figure 3 (Core vs Tail gap profile)
│
├── figures/                       # Publication figures (300 DPI PNG & vector PDF)
│   ├── Fig1_34tests_300dpi.png
│   ├── Fig1_34tests_300dpi.pdf
│   ├── Fig2_length_PGLS.png
│   ├── Fig2_length_PGLS.pdf
│   ├── Fig3_tail_localisation_v3_final.png
│   └── Fig3_tail_localisation_v3_final.pdf
│
└── tables/                        # Statistical summary tables
    ├── Table2_protein_lengths_tail.csv
    ├── Table_all_34_statistical_tests.csv
    └── Table_OLS_PGLS_summary.csv
```

---

## Quickstart & Reproduction

### 1. Installation

Clone the repository and install the Python dependencies:

```bash
git clone https://github.com/Sudheesh02/chaperone-mrna-thermal-adaptation.git
cd chaperone-mrna-thermal-adaptation
pip install -r requirements.txt
```

*(Alternatively, create a conda environment using `conda env create -f environment.yml && conda activate chaperone-adaptation`)*.

### 2. Run the Full Reproduction Pipeline

To re-run all analytical steps, recalculate all test statistics, regenerate all figures, and automatically assert all findings against the manuscript:

```bash
python run_reproduction_pipeline.py
```

### Expected Output Summary:

```text
==================================================================================
 REPRODUCIBILITY PIPELINE: BACTERIAL CHAPERONE THERMAL ADAPTATION
 Manuscript: bioRxiv preprint (Sudheesh et al.)
==================================================================================

[1/5] RUNNING: Mann-Whitney tests, Cliff's delta & BH-FDR / Bonferroni (01_run_34_tests.py)...
✓ Completed in 3.41s

[2/5] RUNNING: Phylogenetic Generalized Least Squares (PGLS) (02_pgls_regression.py)...
✓ Completed in 11.49s

[3/5] RUNNING: GroEL protein alignment & C-terminal tail localization (03_tail_localization.py)...
✓ Completed in 2.47s

[4/5] RUNNING: Altschul-Erickson dinucleotide shuffle verification (04_dinuc_shuffle.py)...
✓ Completed in 0.04s

[5/5] GENERATING PUBLICATION FIGURES & GRAPHICAL ABSTRACT...
  ✓ plot_figure1.py (5.28s)
  ✓ plot_figure2.py (4.46s)
  ✓ plot_figure3.py (3.41s)
  ✓ plot_graphical_abstract.py (5.09s)

==================================================================================
 VERIFYING COMPUTED VALUES AGAINST MANUSCRIPT CLAIMS:
==================================================================================
 [PASS] 34-test battery: exactly 2 comparisons survive BH FDR 5% (groL length & groL z-55°C).
 [PASS] Continuous Topt scaling & PGLS: groL length decreases with Topt (p < 1e-4, lambda ~ 1.0); dnaK invariant.
 [PASS] Protein-level localization: GroEL full length (544.4 vs 541.2 aa) and tail (21.8 vs 18.9 aa) confirmed.
 [PASS] Publication figures: Figures 1, 2, and 3 successfully verified in figures/ directory.

==================================================================================
 Pipeline complete: all manuscript assertions successfully validated.
==================================================================================
```

---

## Individual Analysis Scripts

### 1. `scripts/01_run_34_tests.py`
Computes two-sided Mann-Whitney $U$ tests, Cliff's $\delta$, Benjamini-Hochberg FDR ($\alpha = 0.05$), and Bonferroni corrections across all 34 pre-specified comparisons:
- 15 thermodynamic descriptors $\times$ 2 genes (*dnaK*, *groL*) folded at 37°C and 55°C.
- 4 composition-controlled structural $z$-scores ($N = 100$ Altschul-Erickson shuffles).

### 2. `scripts/02_pgls_regression.py` & `scripts/02_pgls_regression.R`
Fits ordinary least squares (OLS) and phylogenetic generalized least squares (PGLS) regressions of coding sequence length against $T_{opt}$.
- Pagel's $\lambda$ is estimated by maximum likelihood over $\lambda \in [0, 1]$ using the 16S rRNA maximum-likelihood tree (`sequences/16S_ML_tree.nwk`).
- **Software Provenance:** Primary reported manuscript parameters and Figure 2 curves are derived from R `nlme::gls` (`corPagel`: Pagel's $\lambda = 0.997$, slope $= -0.356\text{ nt/}^\circ\text{C}$, $P = 6.0 \times 10^{-5}$). The companion Python script provides an independent implementation that confirms the same qualitative and statistical result ($\lambda \approx 1.0$, slope $= -0.352\text{ nt/}^\circ\text{C}$, $P = 1.6 \times 10^{-5}$) within standard optimizer convergence tolerance. Both model outputs are compiled in `tables/Table_OLS_PGLS_summary.csv`.

### 3. `scripts/03_tail_localization.py`
Reads the 569-column MAFFT alignment of GroEL protein sequences, computes column-by-column gap frequencies in mesophiles and thermophiles, calculates excess deletion rates, and measures excess gap mass in the disordered C-terminal tail (columns $\ge 540$, corresponding to residues 524–548 in *E. coli* numbering).

### 4. `scripts/04_dinuc_shuffle.py`
Stand-alone implementation of the Altschul & Erickson (1985) dinucleotide-preserving sequence permutation algorithm. Preserves exact mononucleotide, dinucleotide, and nearest-neighbor stacking thermodynamics while randomizing higher-order sequence patterns.

### 5. `scripts/plot_figure1.py`, `plot_figure2.py`, `plot_figure3.py`
Regenerates publication Figures 1, 2, and 3 at 300 DPI PNG and vector PDF formats matching the bioRxiv manuscript.

---

## Mapping Data Files to Manuscript Claims

| Paper Finding / Claim | Primary Data File | Script | Key Reported Metric |
|---|---|---|---|
| **50-Strain Cohort** (25 meso, 25 thermo; 4 phyla) | `data/strain_metadata.csv` | `01_run_34_tests.py` | 50 strains, balanced across 4 phyla |
| **34 Statistical Tests** | `data/*_folding_descriptors.csv`, `data/*_zscores_N100.csv` | `01_run_34_tests.py` | 2 survive BH-FDR (*groL* length, *groL* $z_{55}$) |
| **Residual *groL* $z$-Score** at 55°C | `data/groL_zscores_N100.csv` | `01_run_34_tests.py` | $z = -1.41$ vs $-2.77$, $P = 6.4 \times 10^{-4}$ |
| ***groL* CDS Length Shortening** | `data/groL_folding_descriptors.csv` | `01_run_34_tests.py` | $1636.3$ vs $1626.5$ nt, $P = 5.1 \times 10^{-4}$ |
| **Continuous Scaling with $T_{opt}$** | `data/groL_folding_descriptors.csv` | `02_pgls_regression.py` | Spearman $\rho = -0.589$, $P = 6.8 \times 10^{-6}$ |
| **PGLS Phylogenetic Correction** | `sequences/16S_ML_tree.nwk` | `02_pgls_regression.py` | Slope $= -0.356\text{ nt/}^\circ\text{C}$, $\lambda = 0.997$, $P = 6.0 \times 10^{-5}$ |
| **Tail Gap Mass Concentration** | `data/groL_column_gap_profile.csv` | `03_tail_localization.py` | $3.00 / 4.08 = 73.5\%$ in tail ($\ge \text{col } 540$) |
| **GroEL Protein Length Difference** | `data/protein_length_and_tail.csv` | `03_tail_localization.py` | $544.4$ vs $541.2$ aa ($P = 5.1 \times 10^{-4}$) |
| ***dnaK* Invariant Control** | `data/dnaK_folding_descriptors.csv` | `01_run_34_tests.py`, `02_pgls_regression.py` | CDS length $P = 0.71$, protein length $P = 0.71$ |

---

## Citation

```bibtex
@article{singh2026chaperone,
  title={Chaperone mRNA secondary structure does not strongly track thermal adaptation in bacteria, but groL coding sequence shortens in its disordered C-terminal tail},
  author={Singh, Sudheesh and Sahu, Vaibhav and Sradha, Sushree and Dhiman, Vedikaa},
  journal={bioRxiv},
  year={2026}
}
```

---

## License

This repository is licensed under the [MIT License](LICENSE).  
The datasets and sequence files are made available under public domain / [Creative Commons Attribution 4.0 International (CC-BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
