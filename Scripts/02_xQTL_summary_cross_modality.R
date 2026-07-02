################################
# Cross-modality xQTL mapping visualization
################################



# Load Libraries
library(data.table)
library(tidyverse)
library(eulerr)
library(biomaRt)



## -----------------------------
## Colors & Themes
## -----------------------------
colors_modality <- c( "Expression" = "#ff595e", "Isoforms" = "#1982c4", "Splicing" = "#8ac926", "Alt_TSS" = "#9f4cfd", "Alt_polyA" = "#ff924c", "Stability" = "#ffdd47" )

theme_paper <- function(base_size = 12) {
  theme_classic() +
    theme(
      axis.text = element_text(color = "black", face = "bold", size = base_size),
      axis.title = element_text(face = "bold", size = base_size + 2),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(face = "bold", size = base_size, color = "black")
    )
}



## -----------------------------
## Load xQTL files for each transcriptomic modality + cross-modality
## -----------------------------
files <- c(
  Expression = "path/to/expression.cis_independent_qtl.txt.gz",
  Stability  = "path/to/stability.cis_independent_qtl.txt.gz",
  Isoforms   = "path/to/isoforms.cis_independent_qtl.txt.gz",
  Splicing   = "path/to/splicing.cis_independent_qtl.txt.gz",
  Alt_TSS    = "path/to/alt_TSS.cis_independent_qtl.txt.gz",
  Alt_polyA  = "path/to/alt_polyA.cis_independent_qtl.txt.gz",
  Cross = "path/to/cross_modality.cis_independent_qtl.txt.gz"
)

modalities <- imap(files, function(file, mod) {
  dat <- fread(file)
  dat %>% 
    mutate(
      gene_id = if (mod %in% c("Expression", "Stability")) phenotype_id else group_id, 
      modality = mod)
})



## -----------------------------
## Figure 2A
## Number of xGenes by the number of independent xQTL signals (Separate- vs Cross- Modality Mapping)
## -----------------------------
sep_gene_counts <-
  modalities[names(modalities) != "Cross"] %>%
  bind_rows() %>%
  count(gene_id, name = "n_qtl") %>%
  mutate(method = "Separate")

cross_gene_counts <-
  modalities$Cross %>%
  count(gene_id, name = "n_qtl") %>%
  mutate(method = "Combined")

gene_counts <-
  bind_rows(sep_gene_counts, cross_gene_counts) %>%
  mutate(
    bin = factor(if_else(n_qtl >= 8, "8+", as.character(n_qtl)), levels = c(as.character(1:7), "8+")),
  ) %>%
  count(method, bin, name = "n_genes")

# Plot
plot_2a <- 
  ggplot(gene_counts, 
         aes(x = bin, y = n_genes, fill = method)) +
  geom_col(position = "dodge", color = "black") + 
  scale_fill_manual(values = c(Separate = "#E69F00", Combined = "#2c7fb8" )) +
  labs(x = "Total xQTLs per gene", y = "xGenes", fill = "Modality mapping\nmethod") +
  theme_paper(13) +
  theme(
    legend.position = c(0.75, 0.80),
    legend.background = element_rect(fill = "transparent", color = NA),
  )



## -----------------------------
## Figure 2B, 2C
## Venn diagram: Separate vs Cross-modality 
## -----------------------------
modalities$Cross <- modalities$Cross %>%
  mutate(
    modality = recode(
      sub(":.*", "", phenotype_id), 
      expression = "Expression", 
      stability  = "Stability", 
      isoforms   = "Isoforms", 
      splicing   = "Splicing", 
      alt_TSS    = "Alt_TSS", 
      alt_polyA  = "Alt_polyA")
  )

sep_qtls <-
  modalities[names(modalities) != "Cross"] %>%
  bind_rows() %>%
  select(gene_id, modality)

cross_qtls <-
  modalities$Cross %>%
  select(gene_id, modality)

get_venn_objects <- function(sep_genes, cross_genes){
  counts <- c(
    length(setdiff(sep_genes, cross_genes)),
    length(setdiff(cross_genes, sep_genes)),
    length(intersect(sep_genes, cross_genes))
  )
  
  fit <- euler(c(
    "separate" = counts[1], 
    "combined" = counts[2], 
    "separate&combined" = counts[3]))
  
  q_labels <- sprintf("%s\n(%.0f%%)", format(counts, big.mark = ","), counts /sum(counts) * 100)
  
  list(fit = fit, labels = q_labels)
}


venn_plot <- function(sep_df, cross_df, mode = c("all", "per-modality"), mod_name = NULL, colors_modality){
  
  # ----------
  # Summed over all modalities
  # ----------
  if (mode == "all"){
    sep_genes <- unique(sep_df$gene_id)
    cross_genes <- unique(cross_df$gene_id)
    res <- get_venn_objects(sep_genes, cross_genes)
    
    plot(res$fit, 
      quantities = list(labels = res$labels, cex = 1.5, font = 2),
      labels = FALSE,
      legend = list(side = "bottom", cex = 1.2, font = 2, nrow = 1, ncol = 2),
      fills = list(fill = c("#E5E5E5", "#999999"), alpha = 0.8),
      main = list(label = "xGenes per mapping method", cex = 1.5, font = 2)
    ) 
  }
  
  # ----------
  # For each modality
  # ----------
  else{
    sep_genes  <- sep_df %>% filter(modality == mod_name) %>% pull(gene_id) %>% unique()
    cross_genes <- cross_df %>% filter(modality == mod_name) %>% pull(gene_id) %>% unique() 
    res <- get_venn_objects(sep_genes, cross_genes)
    
    plot(res$fit,
         quantities = list(labels = res$labels, cex = 1, font = 2),
         labels = FALSE,                                
         legend = list(side = "bottom", cex = 1, font =2, nrow = 1, ncol = 2),       
         fill = c(
           scales::alpha(colors_modality[mod_name], 0.3), 
           scales::alpha(colors_modality[mod_name], 0.6), 
           colors_modality[mod_name]                      
         ),
         main = list(label = mod_name, cex = 1, font = 2),
    )  
  }
}

# Example: Figure 2B
venn_plot(sep_qtls, cross_qtls, mode = "all", colors_modality = colors_modality)
# Example: Figure 2C - Expression
venn_plot(sep_qtls, cross_qtls, mode = "per-modality", mod_name = "Expression", colors_modality = colors_modality)



## -----------------------------
## Positional distributions of xQTLs 
## -----------------------------

# 1. Prepare data

### Extract variant position
extract_pos <- function(df) {
  df[, variant_pos := as.numeric(tstrsplit(variant_id, ":", keep = 2)[[1]])]
  df
}
modalities <- lapply(modalities, extract_pos)

### Retrieve gene start and end position
### Connect to Ensembl to retrieve gene coordinates
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

get_gene_pos <- function(qtl_df, mart){
  as.data.table(getBM(
    attributes = c('ensembl_gene_id', 'start_position', 'end_position', 'strand'),
    filters = 'ensembl_gene_id',
    values = unique(qtl_df$gene_id),
    mart = mart
  ))
}

### Compute normalized position of variants to its gene length
get_norm_pos <- function(qtl_df, gene_info){
  
  df <- merge(qtl_df, gene_info, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
  
  df[, gene_length := end_position - start_position]
  
  df[, dist_from_tss := fcase(
    strand == 1, variant_pos - start_position, 
    strand == -1, end_position - variant_pos)]
  
  df[, norm_pos := fcase(
    dist_from_tss < 0, dist_from_tss / gene_length,
    dist_from_tss > gene_length, 1 + ((dist_from_tss - gene_length) / gene_length), gene_length == 0, 0,
    dist_from_tss >= 0 & dist_from_tss <= gene_length, dist_from_tss / gene_length
  )]
  
  df[norm_pos >= -1 & norm_pos <= 2]
}

# 2. Histogram generator
make_hist_plot <- function(df, title, fill_color){
  
  ggplot(df, aes(x = norm_pos)) +
    geom_histogram(bins = 100, fill = fill_color, color = NA) +
    geom_vline(xintercept = c(0, 1), color = "gray40", linewidth = 0.6) +
    scale_x_continuous(limits = c(-1, 2), breaks = c(-1, 0, 1, 2), labels = c("", "Gene start", "Gene end", ""), expand = c(0, 0)) +
    labs(title = title) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(color = "black", size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0)
    )
}


# 3. remain/remove plot pair generator
remain_remove_plot <- function(plot_data, modality){
  
  cross_list <- split(modalities$Cross$variant_id, modalities$Cross$modality)
  variant_in_cross <- cross_list[[modality]]
  
  remain_df <- plot_data[variant_id %in% variant_in_cross]
  remove_df <- plot_data[!variant_id %in% variant_in_cross]
  
  remain_plot <- make_hist_plot(remain_df, paste(modality, "(Remained)"), colors_modality[[modality]])
  remove_plot <- make_hist_plot(remove_df, paste(modality, "(Removed)"), colors_modality[[modality]])
  
  # use a common y-axis scale for remained and removed plots
  max_y <- max(ggplot_build(remain_plot)$data[[1]]$count, ggplot_build(remove_plot)$data[[1]]$count, na.rm = TRUE)
  
  remain_plot <- remain_plot +
    scale_y_continuous(
      limits = c(0, max_y * 1.01),
      breaks = pretty(c(0, max_y), n = 2),
      expand = expansion(mult = c(0, 0.05))
    )
  
  remove_plot <- remove_plot +
    scale_y_continuous(
      limits = c(0, max_y * 1.01),
      breaks = pretty(c(0, max_y), n = 2),
      expand = expansion(mult = c(0, 0.05))
    )
  
  list(remain_plot, remove_plot)
}

qtl_position_plot <- function(modality, mart, mode = c("separate", "remain-remove")){
  mode <- match.arg(mode)
  gene_info <- get_gene_pos(modalities[[modality]], mart)
  plot_data <- get_norm_pos(modalities[[modality]], gene_info)
  
  # Positional distributions of xQTLs from separate modality mapping
  if (mode == "separate"){
    make_hist_plot(plot_data, modality, colors_modality[[modality]])
  }
  
  # Positional distributions of xQTLs that were identified by separate-modality mapping and retained, versus the redundant ones removed by cross-modality mapping.
  else{
    remain_remove_plot(plot_data, modality)
  }
}

# Example
qtl_position_plot("Expression", mart, mode = "separate")
qtl_position_plot("Isoforms", mart, mode = "remain-remove")
