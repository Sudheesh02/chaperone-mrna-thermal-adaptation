#!/usr/bin/env Rscript
# 02_pgls_regression.R
# --------------------
# PGLS of groL CDS length vs Topt using R ape + nlme (corPagel).
# Run from the repository root:
#     Rscript scripts/02_pgls_regression.R
# or from scripts/:
#     Rscript 02_pgls_regression.R

suppressPackageStartupMessages({
  if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape", repos="https://cloud.r-project.org")
  if (!requireNamespace("nlme", quietly = TRUE)) install.packages("nlme", repos="https://cloud.r-project.org")
  library(ape)
  library(nlme)
})

# Find repository root
args <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_name) > 0) {
  script_dir <- dirname(normalizePath(script_name))
  base_dir <- dirname(script_dir)
} else {
  base_dir <- getwd()
}

pheno_file <- file.path(base_dir, "data", "groL_folding_descriptors.csv")
tree_file  <- file.path(base_dir, "sequences", "16S_ML_tree.nwk")

pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
tree  <- read.tree(tree_file)

pheno$tip <- gsub("[^A-Za-z0-9]+", "_", pheno$species)
pheno$tip <- gsub("_+", "_", pheno$tip)
common <- intersect(tree$tip.label, pheno$tip)

if (length(common) < nrow(pheno)) {
  for (i in seq_len(nrow(pheno))) {
    if (!(pheno$tip[i] %in% tree$tip.label)) {
      hit <- grep(pheno$tip[i], tree$tip.label, value = TRUE, fixed = TRUE)
      if (length(hit) == 1) pheno$tip[i] <- hit
    }
  }
  common <- intersect(tree$tip.label, pheno$tip)
}

tree2  <- keep.tip(tree, common)
pheno2 <- pheno[pheno$tip %in% common, ]
rownames(pheno2) <- pheno2$tip
pheno2 <- pheno2[tree2$tip.label, ]

cat("Fitting PGLS with Pagel's lambda for groL CDS length vs Topt (n = 50)...\n")
lambda_fit <- gls(cds_length_nt ~ topt_c,
                  data = pheno2,
                  correlation = corPagel(1, tree2, form = ~1, fixed = FALSE),
                  method = "ML")
print(summary(lambda_fit))
cat("\nEstimated Pagel's lambda:", as.numeric(lambda_fit$modelStruct$corStruct), "\n\n")

cat("Fitting Ordinary Least Squares (OLS)...\n")
ols <- lm(cds_length_nt ~ topt_c, data = pheno2)
print(summary(ols))
