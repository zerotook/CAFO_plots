# =========================
# ARG Season Comparison Plot (Styled like the example figure)
# Fold column = Fall_mean / Spring_mean (NOT log2)
# Dot size = fold value
# =========================

library(tidyverse)
library(patchwork)

# -------- SETTINGS --------
N_TOP <- 30          # set to 15 if you want top 15, or Inf for all
OUT_W <- 13
OUT_H <- 9

# 1) Read CSV
df <- readr::read_csv(
  "inputs/ARGsOAP_RA.csv",
  show_col_types = FALSE
)

names(df) <- trimws(names(df))
stopifnot("type" %in% names(df))

# 2) Pick raw replicate columns ending in _F or _S
exclude_cols <- c("fall D","FALL SD","spring D","spring SD")
rep_cols <- setdiff(names(df), exclude_cols) %>%
  keep(~ stringr::str_detect(.x, "_[FS]$"))

# OPTIONAL: exclude UI/WW if you don’t want them treated as replicates
# rep_cols <- rep_cols %>% discard(~ .x %in% c("UI_F","WW_F","UI_S","WW_S"))

# 3) Long format + site type
long <- df %>%
  select(type, all_of(rep_cols)) %>%
  pivot_longer(-type, names_to = "sample", values_to = "value") %>%
  mutate(
    season = if_else(str_ends(sample, "_F"), "Fall", "Spring"),
    season = factor(season, levels = c("Fall","Spring")),
    value  = as.numeric(value),
    sample_prefix = str_remove(sample, "_[FS]$"),
    # ---- EDIT THIS IF NEEDED ----
    site_type = case_when(
      str_detect(sample_prefix, "^D") ~ "Dairy",
      str_detect(sample_prefix, "^S") ~ "Swine",
      TRUE ~ "Other"
    )
  ) %>%
  filter(!is.na(value), site_type %in% c("Dairy","Swine"))

# --------- STYLE HELPERS ---------
zebra_layer <- function(type_levels, alpha = 0.18) {
  n <- length(type_levels)
  tibble(
    ymin = seq(1, n, by = 2) - 0.5,
    ymax = seq(1, n, by = 2) + 0.5,
    xmin = -Inf, xmax = Inf
  ) %>%
    geom_rect(
      data = .,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey70",
      alpha = alpha
    )
}

theme_matrix <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "grey85"),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 9, color = "black"),
      plot.title   = element_text(face = "bold", size = 12),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# 4) Main function: build styled 3-column plot for one site type
make_styled_plot <- function(long_site, site_label, out_png) {
  
  # Stats per type
  stats <- long_site %>%
    group_by(type) %>%
    summarise(
      mean_fall   = mean(value[season == "Fall"], na.rm = TRUE),
      mean_spring = mean(value[season == "Spring"], na.rm = TRUE),
      p_value = {
        x <- value[season == "Fall"]; y <- value[season == "Spring"]
        if (length(x) >= 2 && length(y) >= 2) wilcox.test(x, y)$p.value else NA_real_
      },
      .groups = "drop"
    )
  
  # Pseudocount for dividing by (near) zero spring means
  min_nonzero <- long_site %>% filter(value > 0) %>% summarise(m=min(value)) %>% pull(m)
  pc <- ifelse(is.finite(min_nonzero), min_nonzero/2, 1e-6)
  
  stats <- stats %>%
    mutate(
      fold = (mean_spring + pc) / (mean_fall + pc) ,      # <-- Fall / Spring
      mean_abund = (mean_fall + mean_spring)/2,
      sig = case_when(
        is.na(p_value)  ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  # stats <- stats %>%
  #   mutate(
  #     fold_raw = (mean_fall + pc) / (mean_spring + pc),
  #     fold = pmin(fold_raw, 10)   # cap at 10 for plotting
  #   )
  
  # order types by mean abundance
  type_order <- stats %>% arrange(mean_abund) %>% pull(type)
  if (is.finite(N_TOP)) type_order <- head(type_order, N_TOP)
  
  long_site <- long_site %>%
    filter(type %in% type_order) %>%
    mutate(type = factor(type, levels = type_order))
  
  stats <- stats %>%
    filter(type %in% type_order) %>%
    mutate(type = factor(type, levels = type_order))
  
  type_levels <- levels(long_site$type)
  
  # ---- Panel A: Boxplots ----
  # ---- Panel A: Boxplots ----
  p_boxes <- ggplot(long_site, aes(x = value, y = type, fill = season)) +
    zebra_layer(type_levels, alpha = 0.18) +
    geom_boxplot(
      position = position_dodge2(width = 0.72, preserve = "single"),
      width = 0.62,
      outlier.size = 0.9,
      linewidth = 0.25,
      alpha = 0.95
    ) +
    scale_fill_manual(
      values = c("Fall" = "#ffa384", "Spring" = "#b5e5cf"),
      breaks = c("Fall", "Spring"),
      name = NULL
    ) +
    #labs(title = paste0(site_label, "  |  Seasonal distribution"), x = "Relative abundance") +
    theme_matrix() +
    theme(legend.position = "right")
  p_boxes <- p_boxes +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  
  # ---- Panel B: p stars (NO vertical grid lines) ----
  p_sig <- ggplot(stats, aes(x = 1, y = type)) +
    zebra_layer(type_levels, alpha = 0.18) +
    geom_text(aes(label = sig), size = 4, color = "black") +
    coord_cartesian(xlim = c(0.9, 1.1)) +
    labs(title = "p") +
    theme_matrix() +
    theme(
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),          # <- removes vertical lines + all grid in Panel B
      panel.grid.major.x = element_blank()   # <- extra safety
    )
  
  
  # ---- Panel C: Fold dots (NO x-axis). Size + fill gradient; black border; 0 = hollow ----
  stats <- stats %>%
    mutate(
      fold_cap  = pmin(pmax(fold, 0), 10),
      fold_fill = if_else(fold_cap == 0, NA_real_, fold_cap)  # NA -> hollow
    )
  
  p_fold <- ggplot(stats, aes(x = 1, y = type)) +
    zebra_layer(type_levels, alpha = 0.18) +
    geom_point(
      aes(size = fold_cap, fill = fold_fill),
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.95
    ) +
    coord_cartesian(xlim = c(0.9, 1.1)) +
    scale_size_continuous(
      limits = c(0, 10),
      range  = c(1.5, 6),
      breaks = c(0, 2.5, 5, 7.5, 10),
      name   = "fold"
    ) +
    scale_fill_gradient(
      low = "white",            # small values light
      high = "darkred",         # large values dark  (pick any dark color you like)
      limits = c(0, 10),
      breaks = c(0, 2.5, 5, 7.5, 10),
      na.value = "transparent", # fold=0 hollow
      name = "fold"
    ) +
    theme_matrix() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      panel.grid.major.x = element_blank(),   # <- extra safety
      legend.position = "bottom"
    )
  
  
  
  
  final <- p_boxes + p_sig + p_fold +
    plot_layout(widths = c(3.8, 0.55, 1.8))
  
  ggsave(out_png, final, width = OUT_W, height = OUT_H, dpi = 300)
  final
}

# 5) Build Dairy + Swine
p_dairy <- make_styled_plot(
  long_site  = long %>% filter(site_type == "Dairy"),
  site_label = "Dairy",
  out_png    = file.path("examples", "ARG_style_matrix_Dairy_fold.png")
)

p_swine <- make_styled_plot(
  long_site  = long %>% filter(site_type == "Swine"),
  site_label = "Swine",
  out_png    = file.path("examples", "ARG_style_matrix_Swine_fold.png")
)

p_dairy
p_swine
