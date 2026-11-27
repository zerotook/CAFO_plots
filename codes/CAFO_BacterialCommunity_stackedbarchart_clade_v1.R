library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load data
  df <- read_csv("inputs/Class_clade_clean.csv")

# Long format and filter by abundance threshold per sample
long_df <- df %>%
  select(-lineage) %>%
  pivot_longer(cols = -name, names_to = "site", values_to = "abundance") %>%
  filter(abundance > 0.01)

# Optional: make site order consistent
long_df$site <- factor(long_df$site, levels = unique(long_df$site))

# Plot
p <- ggplot(long_df, aes(x = site, y = abundance, fill = name)) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(
    title = "Bacterial Taxonomic Profiles (Species Level)",
    x = "Sampling Site",
    y = "Abundance",
    fill = "Class"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

print(p)
ggsave("examples/taxonomic_profiles_species_gt_0.05.png", p, width = 12, height = 6, dpi = 300)
