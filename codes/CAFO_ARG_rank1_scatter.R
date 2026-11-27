library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# 1. Read the CSV (your file just converted from Excel)
df <- read.csv(
  "inputs/Rank1_scatterplot.csv",
  check.names = FALSE
)

# First column has no header in the CSV, so give it a proper name
names(df)[1] <- "ARG_class"

# 2. Identify sample columns (all except the ARG_class label)
sample_cols <- setdiff(names(df), "ARG_class")

# 3. Split into:
#    - arg_rows: aminoglycoside, sulfonamide, tetracycline
#    - total_row: "total risk 1" (your total rank 1)
arg_rows <- df$ARG_class %in% c("aminoglycoside", "sulfonamide", "tetracycline")
total_row <- df$ARG_class == "total risk 1"

arg_df <- df[arg_rows, ]
total_df <- df[total_row, ]
# for (var_name in c("sulfonamide", "tetracycline", "aminoglycoside")) {
#   sub <- plot_df %>% filter(ARG_class == var_name)
#   cat(var_name, ":\n")
#   print(shapiro.test(sub$ARG_value))
# }
# cat("total_rank1:\n")
# print(shapiro.test(plot_df$total_rank1[1:20]))
# 4. Long format for ARG values
arg_long <- arg_df %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "ARG_value"
  )

# 5. Long format for total rank 1
total_long <- total_df %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "total_rank1"
  ) %>%
  select(sample, total_rank1)

# 6. Join so each row has: sample, ARG_class, ARG_value, total_rank1
plot_df <- arg_long %>%
  left_join(total_long, by = "sample")

# 7. Plot in your requested style
ggplot(plot_df, aes(x = ARG_value, y = total_rank1)) +
  geom_point(aes(color = ARG_class, shape = ARG_class), size = 3, alpha = 0.8) +
  geom_rug(aes(color = ARG_class), sides = "b", alpha = 0.5) +
  geom_smooth(
    aes(color = ARG_class),
    method    = "lm",
    se        = FALSE,
    fullrange = TRUE
  ) +
  scale_color_manual(
    values = c(
      "sulfonamide"    = "#F5A834",
      "tetracycline"   = "#D62728",
      "aminoglycoside" = "#AEC7E8"
    )
  ) +
  # ggpubr::stat_cor(
  #   aes(color = ARG_class),
  #   method      = "spearman",      # <-- switched from pearson to spearman
  #   label.x.npc = "top",
  #   label.y.npc = "left",
  #   size        = 5
  # ) +
  labs(
    x = "Total abundance of selected ARGs \n (cell per gene copies)",
    y = "Total Rank 1 abundance (cell per gene copies)",
    color = "ARG class",
    shape = "ARG class"
  ) +
  theme_pubr(border = TRUE) +
  theme(legend.position = "right",
        legend.direction = "vertical")
