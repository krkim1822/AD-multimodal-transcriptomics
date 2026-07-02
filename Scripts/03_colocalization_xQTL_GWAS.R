################################
## Colocalization Analysis using coloc.susie
## between xQTL data and AD GWAS summary data
################################

# Load packages
library(data.table)
library(coloc)
library(tidyverse)



## -----------------------------
## File Paths and Parameters
## -----------------------------
# Pre-filtered GWAS data: contains only SNPs within the cis-windows or risk loci for computational efficiency.
GWAS_FILE <- "/path/to/filtered_AD_gwas_tsv.gz" 
# Replace gwas_n, qtl_n, and s_value with with your GWAS sample size, case proportion of your GWAS data, and xQTL sample size.
GWAS_N <- gwas_n 
S_VALUE <- s_value 
QTL_N <- qtl_n
# Replace the path to the specific xQTL modality file to test.
QTL_FILE <- "/path/to/cis_independent_qtl.txt.gz" 
# Reference genotype prefix
PLINK_BFILE <- "/path/to/reference_genotype_prefix" 

TEMP_LD_DIR <- "LD_temp"
if (!dir.exists(TEMP_LD_DIR)) dir.create(TEMP_LD_DIR)
OUTPUT_FILE <- "/path/to/coloc_results.csv"



## -----------------------------
## Data Loading
## -----------------------------
gwas <- fread(GWAS_FILE)
qtl <- fread(QTL_FILE)

# Get required columns from GWAS data for coloc.susie 
locus_list <- unique(gwas[, .(Target_SNP, Gene, Chr, Pos, Start, End)])



## -----------------------------
## Run coloc.susie iteratively for each risk locus
## -----------------------------
coloc_results <- list()

for (i in 1:nrow(locus_list)) {
  
  current <- locus_list[i, ]
  message(sprintf("[%d/%d] Processing %s", i, nrow(locus_list), current$Gene))
  
  # Subset GWAS and xQTL data for the locus
  gwas_chunk <- gwas[Target_SNP == current$Target_SNP]
  qtl_chunk <- unique(
    qtl[Chr == current$Chr & Pos >= current$Start & Pos <= current$End,],
    by = "variant_id")
  
  if (nrow(qtl_chunk) == 0) {
    message("  Skipped: No xQTL data for gene ", current$Gene) 
    next
    }
  
  merged_data <- merge(gwas_chunk, qtl_chunk, by = "variant_id", suffixes = c("_gwas", "_xqtl"))
  
  if (nrow(merged_data) == 0) {
    message("  Skipped: No overlapping SNPs between GWAS and QTL")
    next
    }
  
  # Extract variant list for LD matrix 
  variant_list <- merged_data$variant_id
  variant_file <- file.path(TEMP_LD_DIR, sprintf("variants_to_extract_%03d.txt", i))
  ld_prefix <- file.path(TEMP_LD_DIR, sprintf("ld_matrix_%03d",i))
  write.table(variant_list, file = variant_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Generate LD matrix using PLINK
  system(sprintf("plink --bfile %s --extract %s --r --matrix --out %s", PLINK_BFILE, variant_file, ld_prefix))
  
  # Read LD matrix
  ld_file <- paste0(ld_prefix, ".ld")
  if (!file.exists(ld_file)) {
    message("LD matrix generation failed.")
    next
  }
  ld_mat <- as.matrix(read.table(ld_file, header = FALSE))
  rownames(ld_mat) <- colnames(ld_mat) <- variant_list
  
  nan_snps <- unique(c(
    rownames(ld_mat)[rowSums(is.na(ld_mat)) > 0],
    colnames(ld_mat)[colSums(is.na(ld_mat)) > 0]
  ))
  
  if (length(nan_snps) > 0) {
    message("  Removing ", length(nan_snps), " SNPs with NaN LD values")
    keep_snps <- setdiff(variant_list, nan_snps)
    merged_data <- merged_data[variant_id %in% keep_snps]
    ld_mat <- ld_mat[keep_snps, keep_snps, drop = FALSE]
  }
  
  if (nrow(ld_mat) < 2) {
    message("LD matrix too small after NaN removal")
    next
    }
  
# Build SuSiE datasets
  dataset1 <- list(
    snp = merged_data$variant_id,
    beta = merged_data$beta,
    varbeta = merged_data$standard_error^2,
    N = GWAS_N,
    MAF = merged_data$effect_allele_frequency,
    s = S_VALUE,
    type = "cc",
    LD = ld_mat
  )
  
  dataset2 <- list(
    snp = merged_data$variant_id,
    beta = merged_data$slope,
    varbeta = merged_data$slope_se^2,
    N = QTL_N,
    MAF = merged_data$af,
    type = "quant",
    LD = ld_mat
  )
  
  # Run SuSiE 
  S3 <- tryCatch({
    runsusie(dataset1)
  }, error = function(e) {
    message("  Error in SuSiE GWAS:", e$message)
    return(NULL)
  })
  
  S4 <- tryCatch({
    runsusie(dataset2)
  }, error = function(e) {
    message("  Error in SuSiE QTL:", e$message)
    return(NULL)
  })
  
  if (is.null(S3) || is.null(S4)) next
  
  # Run coloc.susie
  susie.res <- tryCatch({
    coloc.susie(S3, S4)
  }, error = function(e) {
    message("  Error in coloc.susie:", e$message)
    return(NULL)
  })
  
  if (is.null(susie.res)) next
  
  res <- as.data.table(susie.res$summary)
  res[, Target_SNP := current$Target_SNP]
  res[, Gene := current$Gene]
  coloc_results[[i]] <- res
}



# Complie colocalization results
final_coloc_results <- rbindlist(coloc_results, fill = TRUE)
fwrite(final_coloc_results, OUTPUT_FILE)