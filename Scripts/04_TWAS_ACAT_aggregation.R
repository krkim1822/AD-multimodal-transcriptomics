################################
## Combine TWAS p-value per gene per modality using ACAT
## for modalities that generate multiple RNA traits per gene.
################################



## -----------------------------
## Prerequisite
## Input data must contain a "Gene_ID" column.
##
## For multi-trait modalities (e.g., isoforms, splicing, alt_TSS, alt_polyA),
## FUSION output typically contains an "ID" column with RNA trait identifiers. 
## (ex. ENSG00000000971__grp_2_upstream_ENST00000695969)
##
## Extract the Ensembl Gene ID from the "ID" column before running this script: 
##
## df$Gene_ID <- gsub("__.*", "", df$ID)
## Optionally, map Ensembl Gene IDs to HGNC symbols 
## (e.g., using biomaRt) to create the "Gene_name" column.
## -----------------------------



# Load packages
library(ACAT)
library(dplyr)
library(purrr)
library(readr)



INPUT_DIR <- "path/to/fusion_result_directory"
OUTPUT_DIR <- "path/to/output/ACAT_result"



## -----------------------------
## Aggregate trait TWAS p-values to gene-level using ACAT
## -----------------------------
run_acat <- function(file_path){
  df <- read_tsv(file_path)
  
  df_acat <- df %>%
    mutate(
      # Adjust boundary p-values for stable Cauchy transformation
      TWAS.P_adj = case_when(TWAS.P <= 0 ~ 1e-16, TWAS.P >= 1 ~ 1-1e-16, TRUE ~ TWAS.P)
    ) %>%
    group_by(Gene_ID) %>% 
    summarise(
      CHR = first(CHR), 
      P0 = first(P0), # TSS (identical per gene)
      Gene_name = first(Gene_name), # If your data does not contain Gene_name (HGNC symbol), you can use first(Gene_ID) instead
      ACAT_p = ACAT(TWAS.P_adj),
      .groups = "drop"
    )
  df_acat
}



## -----------------------------
## Iterate through list of modalities that contain multiple traits per gene 
## -----------------------------
modalities <- c("isoforms", "splicing", "alt_TSS", "alt_polyA")

walk(modalities, function(mod){
  file_path <- file.path(INPUT_DIR, paste0(mod, ".tsv")) # Load FUSION association results 
  result <- run_acat(file_path)
  output_file <- file.path(OUTPUT_DIR, paste0(mod, "_acat.tsv"))
  write_tsv(result, output_file)
  })