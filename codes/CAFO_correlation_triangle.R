library(readr)
library(tidyverse)

# 1. Read data -------------------------------------------------------------
df_raw <- read_csv("inputs/Cross_method.csv",
                   show_col_types = FALSE)

# Keep only numeric columns, in original order
df <- df_raw %>%
  select(where(is.numeric))

# 2. Spearman correlation matrix ------------------------------------------
cor_mat <- cor(df, use = "pairwise.complete.obs", method = "spearman")

# 2b. P-value matrix (Spearman cor.test) ----------------------------------
p_mat <- matrix(NA_real_, nrow = ncol(df), ncol = ncol(df))
colnames(p_mat) <- rownames(p_mat) <- colnames(df)

for (i in seq_len(ncol(df))) {
  for (j in seq_len(ncol(df))) {
    test <- suppressWarnings(
      cor.test(df[[i]], df[[j]], method = "spearman", exact = FALSE)
    )
    p_mat[i, j] <- test$p.value
  }
}

# 3. Build upper / lower matrices -----------------------------------------
low <- matrix(NA_real_, nrow = nrow(cor_mat), ncol = ncol(cor_mat))
up  <- matrix(NA_real_, nrow = nrow(cor_mat), ncol = ncol(cor_mat))

up[upper.tri(up)]   <- cor_mat[upper.tri(cor_mat)]
low[lower.tri(low)] <- cor_mat[lower.tri(cor_mat)]

# Same for p-values (upper triangle only)
up_p <- matrix(NA_real_, nrow = nrow(p_mat), ncol = ncol(p_mat))
up_p[upper.tri(up_p)] <- p_mat[upper.tri(p_mat)]

# 4. Pivot to long format for plotting ------------------------------------
lower_dat <- low %>%
  as.data.frame() %>%
  `colnames<-`(colnames(df)) %>%
  mutate(xvar = colnames(df)) %>%
  pivot_longer(cols = -xvar, names_to = "yvar", values_to = "value")

upper_dat <- up %>%
  as.data.frame() %>%
  `colnames<-`(colnames(df)) %>%
  mutate(xvar = colnames(df)) %>%
  pivot_longer(cols = -xvar, names_to = "yvar", values_to = "value")

upper_p_dat <- up_p %>%
  as.data.frame() %>%
  `colnames<-`(colnames(df)) %>%
  mutate(xvar = colnames(df)) %>%
  pivot_longer(cols = -xvar, names_to = "yvar", values_to = "p")

# Join p-values onto upper_dat
upper_dat <- upper_dat %>%
  left_join(upper_p_dat, by = c("xvar", "yvar"))

# Keep factor levels in original order
var_levels <- colnames(df)

lower_dat <- lower_dat %>%
  mutate(
    xvar = factor(xvar, levels = var_levels),
    yvar = factor(yvar, levels = var_levels)
  )

upper_dat <- upper_dat %>%
  mutate(
    xvar = factor(xvar, levels = var_levels),
    yvar = factor(yvar, levels = var_levels)
  )

# 5. Plot: triangle style with r, p-stars, and |r| size legend ------------

ggplot(lower_dat, aes(x = xvar, y = yvar)) +
  geom_tile(fill = NA, color = "grey80") +
  # size mapped to |r| so legend can have multiple positive breaks
  geom_point(aes(fill = value, size = abs(value)), shape = 22) +
  # correlation text (upper triangle)
  geom_text(
    data = upper_dat,
    aes(color = value, label = round(value, 2)),
    size = 3
  ) +
  # p-value stars only for significant pairs (p < 0.05)
  geom_text(
    data = upper_dat %>% filter(!is.na(p), p < 0.05),
    aes(
      label = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ ""
      )
    ),
    color = "black",
    size = 4,
    vjust = 0
  ) +
  # size legend for |r|
  scale_size_continuous(
    name   = "|Spearman r|",
    range  = c(1.5, 6),                     # min/max square size
    breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)     # multiple squares in legend
  ) +
  labs(x = NULL, y = NULL) +
  # color legend for signed r
  scale_fill_gradient2(
    low = "#c54132", mid = "white", high = "#5fa2e8",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Spearman r"
  ) +
  scale_color_gradient2(
    low = "#c54132", mid = "white", high = "#5fa2e8",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Spearman r"
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
