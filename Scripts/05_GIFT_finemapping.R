################################
# GIFT Fine-mapping Analysis and Visualization 
################################
library(data.table)
library(dplyr)
library(GIFT)
library(ggrepel)
library(ggplot2)

## -----------------------------
## Load directories and inputs for region
## -----------------------------
GWAS_FILE <- "path/to/GWAS/summary/statistics"
GWAS_LD_FILE <- "path/to/reference/LD/data"
SNP_LIST <- snp_list # Vector of SNPs in the target region
WEIGHT_LIST <- weight_list # List of FUSION weight 
GENE_NAMES <- gene_names # Vector of gene names in the target region
GWAS_N <- GWAS_sample_size
TWAS_RESULTS <- TWAS_result # TWAS result data with region column  



## -----------------------------
## Preprocess input data
## -----------------------------
# For function details and input requirements, see source code: 
# https://github.com/yuanzhongshang/GIFT/blob/main/R/pre_process.R

# 1. GWAS summary data pre-process
# Align GWAS summary data with the reference LD matrix and SNPs within the region. 
convert <- pre_process_twostage(
    GWASfile   = GWAS_FILE,
    GWASLDfile = GWAS_LD_FILE,
    snplist    = SNP_LIST
  )



# 2. TWAS (FUSION) weights pre-process
# Converts a list of gene weights into a block-diagonal weight matrix required by GIFT.
betax_mat <- weightconvert(WEIGHT_LIST)



## -----------------------------
## Run GIFT 
## -----------------------------
# For model details, arguments, and default settings, see source code:
# https://github.com/yuanzhongshang/GIFT/blob/main/R/GIFT_two_stage_summ.R

gift_result <- GIFT_two_stage_summ(
  betax = betax_mat,
  betay = convert$beta,
  se_betay = convert$se,
  Sigma = convert$LDmatrix,
  n =  GWAS_N,
  gene = GENE_NAMES, 
  in_sample_LD = FALSE)



## -----------------------------
## Visualization 
## Generate a locus-zoom plot showing GWAS, TWAS and GIFT results for a target region
## -----------------------------
plot_gift_locus <- function(target_region, gwas_df, twas_df, gift_result){
  
  twas_sub <- twas_df %>% filter(Region_group == target_region)
  region_data <- gift_result %>% inner_join(twas_sub, by = c("gene" = "Gene_ID"))
  
  chr_region <- region_data$Chr[1]
  start_region <- min(region_data$Region_start)
  end_region <- max(region_data$Region_end)
  gwas_sub <- gwas_df %>% filter(CHR == chr_region & BP >= start_region & BP <= end_region)
  
  # Region-specific Bonferroni threshold
  # (number of TWAS genes within the target region)
  n_gene_region <- nrow(region_data) 
  bonf_thresh <- -log10(0.05 / n_gene_region)
  gift_significant_genes <- region_data %>% filter(p < 0.05/n_gene_region) 
  
  p <- ggplot()+
    
    # GWAS (gray)
    geom_point(data = gwas_sub, aes(x = BP/1e6, y = -log10(P), color = "GWAS"), alpha = 0.5, size = 1) +
    
    # TWAS (blue squares)
    geom_point(data = region_data, aes(x = (Start+End)/2/1e6, y = -log10(P_value),  color = "TWAS"), shape = 15, size = 2.0) +
    
    # GIFT (red diamond)
    geom_point(data = region_data, aes(x = (Start+End)/2/1e6, y = -log10(p),  color = "GIFT"), shape = 18, size = 3.0) +
    
    # GIFT significant Bonferroni line
    geom_hline(yintercept = bonf_thresh, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_label_repel(data = filter(gift_significant_genes, !grepl("^ENSG", external_gene_name)),
                     aes(x = (Start+End)/2/1e6, y = -log10(p), label = external_gene_name),  
                     fill = "white", color = "black", size = 2.5, label.size = 0.3, fontface = "bold",
                     box.padding = 0.3, point.padding = 0.4, segment.color = "black",segment.size = 0.5, min.segment.length = 0)+
    scale_color_manual(
      name = NULL,
      values = c("GWAS" = "grey70",
                 "TWAS" = "dodgerblue4",
                 "GIFT" = "firebrick3"),
      breaks = c("GWAS", "TWAS", "GIFT")
    ) +
    labs(x = paste0("Position on Chr", chr_region, "(Mb)" ),
         y = "-log10(P-value)") +
    theme_classic() + 
    theme(
      legend.position = "bottom", 
      legend.text = element_text(size = 11, face = "bold"), 
      legend.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(face = "bold") 
    )
  
  p
}



# Example 
# TWAS_RESULTS <- expression_twas
# Expression_region_1 <- plot_gift_locus(
#   target_region = "Region_1", 
#   gwas_df = AD_gwas, 
#   twas_df = TWAS_RESUTS, 
#   gift_result = gift_result)
# print(Expression_region_1)