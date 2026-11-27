#install.packages(c("readr", "dplyr", "tidyr", "RColorBrewer", "ggplot2", "stringr", "ggbreak"))

library(readr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(ggbreak)

file <- "inputs/combined_normalized_cell.subtype_circle - 10_11.csv"

# ---- 1. Read and Process Data -------------------------------------------

cat("Reading file...\n")

df <- read_csv(file, show_col_types = FALSE) %>%
  rename(subtype_full = 1)

sample_cols <- setdiff(names(df), "subtype_full")

dat <- df %>%
  mutate(
    type = sub("__.*$", "", subtype_full),
    subtype = sub("^[^_]+__", "", subtype_full),
    value = rowSums(across(all_of(sample_cols)), na.rm = TRUE)
  ) %>%
  filter(is.finite(value), value > 0)

# ---- FILTER: Exclude Biocide/Defensin/Ionophore --------------------------

dat <- dat %>%
  filter(!type %in% c("biocide", "defensin", "ionophore"))

cat("Found", nrow(dat), "subtypes across", length(unique(dat$type)), "types\n")

# ---- 2. Aggregate subtype and class totals -------------------------------

sub_agg <- dat %>%
  group_by(type, subtype) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  arrange(type, desc(value))

type_tot <- sub_agg %>%
  group_by(type) %>%
  summarise(total = sum(value), .groups = "drop") %>%
  arrange(desc(total))

# ---- 3. Keep top antibiotic classes --------------------------------------

top_types <- 15

n_rows <- nrow(type_tot)

type_tot <- type_tot %>%
  head(min(top_types, n_rows))

sub_agg <- sub_agg %>%
  filter(type %in% type_tot$type)

types <- type_tot$type

grand_total <- sum(type_tot$total)

cat("Grand Total Reads:", grand_total, "\n")

# ---- 4. Custom color scheme ----------------------------------------------

custom_colors <- c(
  "beta_lactam"     = "#FF9896",
  "MLS"             = "#9467BD",
  "macrolide"       = "#9467BD",
  "defensin"        = "#2CA02C",
  "polymyxin"       = "#FF7F0E",
  "rifamycin"       = "#8C564B",
  "fosfomycin"      = "#E377C2",
  "tetracycline"    = "#D62728",
  "trimethoprim"    = "#BCBD22",
  "vancomycin"      = "#17BECF",
  "aminoglycoside"  = "#AEC7E8",
  "phenicol"        = "#98DF8A",
  "chloramphenicol" = "#98DF8A",
  "sulfonamide"     = "#F5A834",
  "quinolone"       = "#C5B0D5",
  "multidrug"       = "#EBE9E8",
  "bacitracin"      = "#8DD3C7",
  "mupirocin"       = "#F7B6D2",
  "novobiocin"      = "#DBDB8D"
)

type_col <- custom_colors[types]

missing_idx <- is.na(type_col)

if (any(missing_idx)) {
  missing_types <- types[missing_idx]
  cat("Warning: No custom color found for:", paste(missing_types, collapse = ", "), "\n")
  
  fallback_pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(missing_types))
  type_col[missing_idx] <- fallback_pal
}

names(type_col) <- types

# ---- 5. Prepare data for stacked bar plot --------------------------------

plot_dat <- sub_agg %>%
  left_join(type_tot, by = "type") %>%
  group_by(type) %>%
  arrange(desc(value), .by_group = TRUE) %>%
  mutate(
    type = factor(type, levels = types),
    subtype_rank = rank(-value, ties.method = "first"),
    within_class_frac = value / total,
    global_frac = value / grand_total,
    fill_id = paste(type, subtype, sep = "___")
  ) %>%
  ungroup()

# ---- 6. Label dominant subtypes only -------------------------------------

# Hard rule:
# Do not label subtype if it contributes <10% within its antibiotic class

label_within_threshold <- 0.10   # 10%

plot_dat <- plot_dat %>%
  group_by(type) %>%
  mutate(
    label = if_else(
      within_class_frac >= label_within_threshold,
      subtype,
      NA_character_
    )
  ) %>%
  ungroup()

# ---- 7. Generate subtype colors from class colors ------------------------

fill_cols <- c()

for (t in types) {
  
  d_t <- plot_dat %>%
    filter(as.character(type) == t) %>%
    arrange(desc(value))
  
  n_sub <- nrow(d_t)
  base_col <- type_col[[t]]
  
  if (is.na(base_col)) {
    base_col <- brewer.pal(8, "Set2")[1]
  }
  
  if (n_sub > 1) {
    cols_t <- colorRampPalette(c(base_col, "white"))(n_sub + 2)[1:n_sub]
  } else {
    cols_t <- base_col
  }
  
  names(cols_t) <- d_t$fill_id
  fill_cols <- c(fill_cols, cols_t)
}

plot_dat <- plot_dat %>%
  arrange(type, desc(value)) %>%
  mutate(fill_id = factor(fill_id, levels = unique(fill_id)))

# ---- 8. Define y-axis breaks ---------------------------------------------
# Based on your data:
# most classes <0.1
# polymyxin ~0.2
# multidrug ~0.4
# bacitracin ~1.0

break1 <- c(0.10, 0.17)  # separates small classes from polymyxin
break2 <- c(0.25, 0.35)  # shortens the gap below multidrug
break3 <- c(0.45, 0.85)  # shortens the bacitracin bar

cat("Y-axis breaks:\n")
cat("Break 1:", break1[1], "to", break1[2], "\n")
cat("Break 2:", break2[1], "to", break2[2], "\n")
cat("Break 3:", break3[1], "to", break3[2], "\n")

# ---- 9. Draw stacked bar plot with multiple broken y-axis -----------------

p <- ggplot(plot_dat, aes(x = type, y = value, fill = fill_id)) +
  
  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.25
  ) +
  
  # geom_text(
  #   aes(label = label),
  #   position = position_stack(vjust = 0.5),
  #   size = 3,
  #   lineheight = 0.85,
  #   color = "black",
  #   na.rm = TRUE
  # ) +
  
  scale_fill_manual(values = fill_cols, guide = "none") +
  
  labs(
    x = "Antibiotic class",
    y = "Total ARG abundance"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  
  # Multiple y-axis breaks
  scale_y_continuous(
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_break(
    break1,
    scales = 0.35,
    space = 0.04
  ) +
  scale_y_break(
    break2,
    scales = 0.35,
    space = 0.04
  ) +
  scale_y_break(
    break3,
    scales = 0.30,
    space = 0.04
  )

# Optional: add visual break marks on the affected columns
break_marks <- tibble(
  type = factor(
    c("polymyxin", "multidrug", "bacitracin"),
    levels = types
  ),
  y = c(0.105, 0.255, 0.455),
  label = c("//", "//", "//")
) %>%
  filter(as.character(type) %in% types)

p <- p +
  geom_text(
    data = break_marks,
    aes(x = type, y = y, label = label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "black"
  )

print(p)
# ---- 10. Save figure -----------------------------------------------------

ggsave(
  "examples/ARG_class_subtype_stacked_bar_absolute_broken_axis.png",
  p,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "white"
)

# ggbreak's PDF device can emit an extra blank page. Keep the 600-dpi plot
# above as the source of a single-page PDF so both example formats are clean.
figure_png <- "examples/ARG_class_subtype_stacked_bar_absolute_broken_axis.png"
figure_pdf <- "examples/ARG_class_subtype_stacked_bar_absolute_broken_axis.pdf"
figure_raster <- png::readPNG(figure_png, native = TRUE)
grDevices::pdf(figure_pdf, width = 12, height = 7, bg = "white", useDingbats = FALSE)
grid::grid.newpage()
grid::grid.raster(
  figure_raster,
  width = grid::unit(1, "npc"),
  height = grid::unit(1, "npc"),
  interpolate = FALSE
)
grDevices::dev.off()
