################################
# Generates the combined covariate matrix for xQTL mapping and TWAS models
# by merging Genotype PCs, RNA trait PCs, and additional sample covariates. 
# Adapted from the default Pantry implementation:
# https://github.com/PejLab/Pantry/blob/main/pheast/scripts/covariates.R
################################


## -----------------------------
## Load packages 
## -----------------------------
suppressPackageStartupMessages(library(snpStats))
library(impute)
library(readr)
library(dplyr)






## -----------------------------
## load_geno: Reads PLINK files 
load_geno <- function(filename) {
    geno <- read.plink(filename)
    # Convert to 0,1,2 coding
    geno_mat <- as(geno$genotypes, "numeric")
    t(geno_mat)
}



## -----------------------------
## get_PCs: Calculates Principal Components
get_PCs <- function(df, n_pcs) {
    if (sum(is.na(df)) > 0) {
        df <- impute.knn(df)$data # Expects samples as columns
    }
    df <- df[apply(df, 1, var) != 0, ]
    pca <- prcomp(t(df), center = TRUE, scale = TRUE)
    n_pcs <- min(n_pcs, ncol(pca$x) - 1)
    pcs <- round(pca$x[, 1:n_pcs, drop = FALSE], 6)
    pcs_df <- data.frame(ID = colnames(pcs), t(pcs))
    colnames(pcs_df) <- c("ID", rownames(pcs))  # Column names starting with digits get 'fixed' and must be changed back
    pcs_df
}



## -----------------------------
## Main Execution
args <- commandArgs(trailingOnly = TRUE)
GENO_PREFIX <- args[1]
BED_FILE <- args[2]
N_GENO_PCS <- as.integer(args[3])
N_PHENO_PCS <- as.integer(args[4])
COV_FILE <- args[5] # Custom external covariate TSV (e.g., input/sample_covar.tsv)
OUT_FILE <- args[6]



## -----------------------------
## Process RNA traits PCs
pheno <- read.delim(BED_FILE, check.names = FALSE, row.names = 4)[, -(1:3)]
if (ncol(pheno) < 2) stop("Computing covariate PCs requires more than 1 sample.")
pheno_pcs <- get_PCs(pheno, N_PHENO_PCS)
pheno_pcs$ID <- paste("pheno", pheno_pcs$ID, sep = "_")



## -----------------------------
## Process Genotype PCs
geno <- load_geno(GENO_PREFIX)
geno <- geno[, colnames(pheno)]
geno_pcs <- get_PCs(geno, N_GENO_PCS)
geno_pcs$ID <- paste("geno", geno_pcs$ID, sep = "_")

stopifnot(identical(colnames(geno_pcs), colnames(pheno_pcs)))
covars <- rbind(geno_pcs, pheno_pcs)



## -----------------------------
## Process additional covariates

# Load additional covariates and format rownames
extra_cov <- read_tsv(COV_FILE, col_types = cols()) %>% rename(ID = 1) %>% as.data.frame()
rownames(extra_cov) <- extra_cov$ID
extra_cov$ID <- NULL

# Ensure categorical variables are treated as factors.

cat_cols <- sapply(extra_cov, function(x) is.character(x) || is.factor(x))

if (any(cat_cols)){
  extra_cov[cat_cols] <- lapply(extra_cov[cat_cols], as.factor)
  
  # Perform one-hot encoding 
  formula_str <- paste("~", paste(names(extra_cov)[cat_cols], collapse = " + "))
  encoded_cats <- model.matrix(as.formula(formula_str), data = extra_cov)[, -1, drop = FALSE] 
  
  # Combine numeric variables and encoded categorical variables
  num_vars <- extra_cov[, !cat_cols, drop = FALSE]
  final_df <- cbind(num_vars, encoded_cats)
} else{
  final_df <- extra_cov
}

# Transpose to match covaraites format for Pheast (Rows: Covaraite, Cols: Samples)
final_df_t <- as.data.frame(t(final_df), row.names = F)
final_df_t <- cbind(ID = colnames(final_df), final_df_t)



## -----------------------------
# Merge all covaraites and save 
stopifnot(identical(colnames(covars), colnames(final_df_t)))
covars_all <- rbind(covars, final_df_t)

# Final covariate matrix 
write.table(covars_all, OUT_FILE, sep="\t", quote=FALSE, row.names=FALSE)
