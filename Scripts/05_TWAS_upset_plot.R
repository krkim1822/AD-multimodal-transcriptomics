################################
# Upset plot of TWAS risk genes across six transcriptomic modalities
################################

# Load packages
library(tidyverse)
library(ggupset)



#-------------------------------
# Input: a named list of character vectors containing TWAS risk genes identified for each modality
#-------------------------------

gene_list <- list(
  Expression = c(...),
  Stability = c(...),
  Isoforms = c(...),
  Splicing = c(...),
  Alt_TSS = c(...),
  Alt_polyA = c(...)
)

# Prepare required data format for upset plot
upset_data <- enframe(gene_list, name = "set", value = "Gene_name") %>%
  unnest(Gene_name) %>% 
  mutate(value = 1) %>%
  pivot_wider(names_from = set, values_from = value, values_fill = 0) %>% 
  pivot_longer(-Gene_name, names_to = "set", values_to = "present") %>%
  filter(present == 1) %>%
  group_by(Gene_name) %>%
  summarise(mod_sets = list(set), .groups = "drop")

# Upset plot
upset_plot <- function(
    plot_data, 
    order_by = "freq", 
    title = NULL, 
    xlab = "Modality group", 
    ylab = "Number of independent TWAS risk genes",
    fill_color = "blue4", 
    base_size = 16){
  ggplot(plot_data, aes(x = mod_sets)) +
  geom_bar(fill = fill_color, width = 0.8) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, fontface = "bold", size = 5) +
  scale_x_upset(order_by = order_by) +
  labs(x = xlab, y = ylab, title = title) +
  theme_bw() +
  theme(
    text = element_text(size = base_size, face = "bold"),
    axis.text = element_text(face = "bold", color = "black", size = base_size),
    axis.title = element_text(face = "bold", size = base_size),
    axis.text.x = element_blank(),
    legend.position = "none"
  ) +
  theme_combmatrix(
    combmatrix.label.text = element_text(face = "bold", color = "black", size = base_size)
  )
}

# Generate plot 
upset_plot(upset_data)
