# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

# 1. Read the data
# Ensure the file 'rank_ARG_bubble.csv' is in your working directory
df <- read.csv("inputs/rank_ARG_bubble.csv", check.names = FALSE)
names(df)[1] <- "sites"  # Normalize the UTF-8 BOM on the first header when present.

# 2. Transform to Long Format
df_long <- df %>%
  pivot_longer(
    cols = -sites,
    names_to = "ARG_Type",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0)

# 3. Create Grouping & Sorting Logic (Fall -> Spring)
df_long <- df_long %>%
  mutate(
    Season = case_when(
      str_ends(sites, "_F") ~ "Fall",
      str_ends(sites, "_S") ~ "Spring",
      TRUE ~ "Other"
    )
  )

df_sorted_sites <- df_long %>%
  select(sites, Season) %>%
  distinct() %>%
  arrange(Season, sites) 

df_long$sites <- factor(df_long$sites, levels = df_sorted_sites$sites)

# 4. Clean Labels and Colors
clean_labels <- function(x) { str_remove_all(x, "_F|_S") }
df_long$ARG_Type <- str_replace_all(df_long$ARG_Type, "_", "-")

ARG_COLORS <- c(
  "beta-lactam"       = "#FF9896",
  "MLS" = "#9467BD",
  "defensin"          = "#2CA02C",
  "polymyxin"         = "#FF7F0E",
  "rifamycin"         = "#8C564B",
  "fosfomycin"        = "#E377C2",
  "tetracycline"      = "#D62728",
  "trimethoprim"      = "#BCBD22",
  "vancomycin"        = "#17BECF",
  "aminoglycoside"    = "#AEC7E8",
  "chloramphenicol"   = "#98DF8A",
  "sulfonamide"       = "#F5A834",
  "quinolone"         = "#C5B0D5",
  "multidrug"         = "#EBE9E8",
  "bacitracin"        = "#8DD3C7",
  "mupirocin"         = "#F7B6D2", 
  "novobiocin"        = "#DBDB8D"
)

# 5. Generate the Plot (With Hollow Legend)
ggplot(df_long, aes(x = sites, y = ARG_Type)) +
  geom_point(aes(size = Abundance, color = ARG_Type), alpha = 0.7) +
  scale_color_manual(values = ARG_COLORS, na.value = "grey80") +
  
  scale_x_discrete(position = "top", labels = clean_labels) +
  scale_size(range = c(2, 12)) +
  
  theme_bw() +
  labs(
    x = "Sampling Site", 
    y = "ARG Type",
    size = "Relative\nAbundance",
    color = "ARG Type"
  ) +
  
  # --- THIS IS THE KEY PART FOR HOLLOW LEGEND DOTS ---
  guides(
    size = guide_legend(
      override.aes = list(
        shape = 21,       # Use shape 21 (circle with border)
        fill = NA,        # No fill (transparent)
        color = "black",  # Black border
        stroke = 1        # Border thickness
      )
    )
  ) +
  # ---------------------------------------------------

theme(
  axis.text.x = element_text(angle = 45, hjust = 0, size = 10, face = "bold"),
  axis.text.y = element_text(size = 11, face = "bold"),
  panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
  legend.position = "right"
)
