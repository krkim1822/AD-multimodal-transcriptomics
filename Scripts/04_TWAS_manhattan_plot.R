################################
# TWAS Manhattan plot 
################################


# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)



## -----------------------------
## Manhattan plotting function 
## -----------------------------
manhattan_plot <- function(
    data, 
    plot_title, 
    p_col = "TWAS.P", 
    snp_col = "Gene_name",
    point_size = 1.3,
    chr_colors = c("blue4", "gray"),
    sig_color = "red",
    line_color = "gray40",
    line_type = "dashed",
    box_padding = 0.3,
    max_overlaps = 10){
  
  # Prepare plot data from TWAS results (ACAT results for multi-trait modalities)
  plot_data <- data %>% 
    mutate(
      CHR = as.numeric(gsub("chr", "", as.character(CHR))), 
      BP = P0, 
      P = .data[[p_col]], 
      LABEL = .data[[snp_col]])
  
  # Get plot position from chromosome positions
  pos_df <- plot_data %>% 
    group_by(CHR) %>% 
    summarise(chr_len = max(BP), .groups = "drop") %>% 
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    left_join(plot_data, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot)
  
  # Define x-axis chromosome labels spacing 
  axisdf <- pos_df %>% 
    group_by(CHR) %>% 
    summarize(center=( max(BPcum) + min(BPcum) ) / 2, .groups = "drop" )
  
  # Bonferroni significance threshold
  bonf_thresh <- 0.05 / nrow(data)
  
  # Draw plot
  man_plot <- ggplot(pos_df, aes(x = BPcum, y = -log10(P))) + 
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = point_size) +
    scale_color_manual(values = rep(chr_colors, 22)) +
    
    # Highlight significant TWAS risk genes in red
    geom_point(data = filter(pos_df, P < bonf_thresh), color = sig_color, size = point_size) +
    
    # Label significant genes with official HGNC symbols (Avoid labeling unmapped Ensembl IDs for better visibility)
    geom_text_repel(data = filter(pos_df, P < bonf_thresh & !grepl("^ENSG", LABEL)), 
                    aes(label = LABEL), size = 3, box.padding = box_padding, max.overlaps = max_overlaps, nudge_y = 1, ylim = c(-log10(bonf_thresh), NA) ) +
    
    # Significance threshold line
    geom_hline(yintercept = -log10(bonf_thresh), color = line_color, linetype = line_type) +
    scale_x_continuous(labels = axisdf$CHR, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0.3)) +  
    labs(title = plot_title, x = "Chromosome", y = "-log10(P)")+
    theme_bw() +
    theme( 
      axis.text = element_text(color = "black", size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.line = element_line(color="black"),
      plot.title = element_text(face = "bold"),
    ) 
  man_plot
}



## -----------------------------
## Usage example
## -----------------------------
# Load FUSION results or ACAT aggregated results
# e.g., expression <- read_tsv("path/to/expression_twas_result.tsv")

# FUSION outputs (Expression, Stability)
manhattan_expression <- manhattan_plot(data = expression, plot_title = "Expression", p_col = "TWAS.P", snp_col = "Gene_name")
print(manhattan_expression)

# ACAT aggregated outputs (Isoforms, Splicing, Alt_TSS, Alt_polyA)
manhattan_isoforms <- manhattan_plot(data = isoforms, plot_title = "Isoform ratios", p_col = "ACAT_p", snp_col = "Gene_name")
print(manhattan_isoforms)