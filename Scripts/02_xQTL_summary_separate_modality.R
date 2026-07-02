################################
# Separate-modality xQTL mapping visualization 
################################



# Load Libraries
library(data.table)
library(tidyverse)
library(forcats)



## -----------------------------
## Colors & Themes
## -----------------------------
colors_qtl_count <- c("1" = "#253494", "2" =  "#2c7fb8", "3" = "#41b6c4", "4" = "#a1dab4", "5+" = "#ffffcc")

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
## Load xQTL files for each transcriptomic modality  
## -----------------------------
files <- c(
  Expression = "path/to/expression.cis_independent_qtl.txt.gz",
  Stability  = "path/to/stability.cis_independent_qtl.txt.gz",
  Isoforms   = "path/to/isoforms.cis_independent_qtl.txt.gz",
  Splicing   = "path/to/splicing.cis_independent_qtl.txt.gz",
  Alt_TSS    = "path/to/alt_TSS.cis_independent_qtl.txt.gz",
  Alt_polyA  = "path/to/alt_polyA.cis_independent_qtl.txt.gz"
)

modalities <- imap(files, function(file, mod) {
  dat <- fread(file)
  dat %>% 
    mutate(
      gene_id = if (mod %in% c("Expression", "Stability")) phenotype_id else group_id)
})


## -----------------------------
## Figure 1A
## Number of unique xGenes per modality
## -----------------------------

# Count the number of xQTLs in each gene across modalities
sep_qtls <- map_df(names(modalities), function(mod){
  modalities[[mod]] %>% 
    count(gene_id, name = "n_qtl") %>% 
    mutate(modality = mod)
})

# Gene counts by number of independent xQTLs
gene_counts <- sep_qtls %>%
  mutate(
    qtl_bin = factor(
      if_else(n_qtl >= 5, "5+", as.character(n_qtl)),
      levels = c("1", "2", "3", "4", "5+")
    )
  ) %>%
  count(modality, qtl_bin, name = "n_genes") %>%
  group_by(modality) %>%
  mutate(
    total_genes = sum(n_genes),
    n_genes_k = n_genes / 1000
  ) %>%
  ungroup() %>%
  mutate(modality = fct_reorder(modality, total_genes, .desc = TRUE))

# Plot 
plot_1a <- 
  ggplot(gene_counts, 
         aes(x = n_genes_k, y = modality, fill = qtl_bin)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = colors_qtl_count, name = "xQTLs") +
  labs(x = "xGenes (x1000)", y = "Modality") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10, face = "bold", color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor.x = element_blank()
  )



## -----------------------------
## Figure 1B
## Number of xGenes by the number of independent xQTL signals, summed over six modalities 
## -----------------------------
qtl_ranks <- lapply(names(modalities), function(mod){
  modalities[[mod]] %>% 
    select(gene_id, rank) %>% 
    mutate(modality = mod) }) %>% 
  bind_rows()

gene_qtl_counts <- qtl_ranks %>%
    group_by(gene_id, modality) %>%
  # Number of conditionally independent xQTLs for each gene
    summarise(n_indep_qtl = max(rank), .groups = "drop") %>%
  mutate(num_group = if_else(n_indep_qtl >= 8, "8+", as.character(n_indep_qtl)))

# Plot 
plot_1b <- 
  gene_qtl_counts %>% 
  count(num_group) %>% 
  ggplot(aes(x = num_group, y=n, fill = n)) +
  geom_col(color = "black") +
  scale_fill_gradient(low = "#a1dab4", high = "#253494") +
  scale_x_discrete(limits = c(1:7, "8+")) +
  labs(
    x = "Number of conditionally independent xQTLs per gene",
    y = "Genes")+
  theme_paper()+
  theme(legend.position = c(0.9, 0.6))



## -----------------------------
## Figure 1C
## Modality proportions of xQTLs of each within-xGene rank
## -----------------------------
rank_summary <- qtl_ranks %>%
  mutate(
    rank_group = factor(
      if_else(rank >= 8, "8+", as.character(rank)), 
      levels = c(as.character(1:7), "8+")
      )) %>%
  count(rank_group, modality, name = "n") %>%
  group_by(rank_group) %>%
  mutate(prop = n / sum(n), total = sum(n)) %>%
  ungroup()

rank_order <- rank_summary %>%
  filter(rank_group == "1") %>%
  arrange(n) %>% 
  pull(modality)

rank_summary <- rank_summary %>% 
  mutate(modality = factor(modality, levels = rank_order))

plot_1c <- 
  ggplot(rank_summary, 
         aes(x = rank_group, y = prop, fill = modality)) +
  geom_col(width = 0.9, color = "black", linewidth = 0.1) +
  geom_text(
    aes(x = rank_group, y = 1.05, label = total),
    data = rank_summary %>% distinct(rank_group, total),
    inherit.aes = FALSE, fontface = "bold", size = 4.5
  ) +
  scale_fill_manual(values = colors_modality) +
  labs(
    x = "Rank of xQTL within gene",
    y = "Proportion of xQTLs",
    fill = "Modality"
  ) +
  guides(fill = guide_legend(reverse = FALSE)) + 
  theme_paper(base_size = 13) +
  theme(legend.position = "right")
