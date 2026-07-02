################################
## Figure 4
## Upset plot of colocalized xGenes across six transcriptomic modalities
################################



# Load required packages
library(data.table)
library(tidyverse)
library(ggupset)


# Requires colocalization result for each modality as *_coloc_results
modalities <- list(Expression = expression_coloc_results, Stability = stability_coloc_results, Isoforms = isoforms_coloc_results, Splicing = splicing_coloc_results, Alt_TSS = alt_TSS_coloc_results, Alt_polyA = Alt_polyA_coloc_results)

# Filter colocalized genes with threshold PP.H4.abf > 0.75 
coloc_threshold <- 0.75
coloc_genes <- map(modalities, function(mod){mod %>% filter(PP.H4.abf > coloc_threshold) %>% pull(Gene) %>% unique()})



## -----------------------------
## Format Data for ggupset
## -----------------------------
all_genes <- unique(unlist(coloc_genes))
plot_data <- tibble(Gene = all_genes) %>%
  mutate(
   # For each gene, find which modalities it belongs to 
    mod = map(Gene, function(current_gene){
      names(coloc_genes)[map_lgl(coloc_genes, function(set) current_gene %in% set)]
    })
    )



## -----------------------------
## Upset Plot
## -----------------------------

# Custom color palette for the intersection degree 
color_mapping <- c("5" = "#4B0055", "4" = "#00608E", "3" = "#009796", "2" = "#C0DE35", "1" = "#FDE333")

plot_3 <- 
  ggplot(plot_data, aes(x = mod)) +
  geom_bar(aes(fill = as.factor(after_stat(count)))) +
  scale_x_upset(order_by = "freq") +
  scale_fill_manual(values = color_mapping, breaks = c("5", "4", "3", "2", "1"), name = "Gene Count") +
  labs(x = "Genes per modality", y = "Number of xGenes") +
  theme_bw() + 
  theme(
    text = element_text(size=16,face = "bold"),
    axis.text = element_text(face = "bold", color = "black", size = 16),
    axis.title = element_text(face = "bold", size = 20),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
    ) +
  theme_combmatrix(
    combmatrix.label.text = element_text(face = "bold", color = "black", size = 16)
    )

print(plot_3)