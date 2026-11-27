library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

df <- read_excel("inputs/qplot.xlsx")

df <- df %>%
  mutate(
    site_id = site_type,
    site_type = ifelse(site_type %in% c("WW", "UI"), "Reference", site_type)
  )

df$season <- factor(df$season, levels = c("fall", "spring"))
df$site_type <- factor(df$site_type, levels = c("Dairy", "Swine & Dairy", "Reference"))
df$site_id <- factor(df$site_id, levels = c("Dairy", "Swine & Dairy", "UI", "WW"))

season_colors <- c("fall" = "#ffa384", "spring" = "#b5e5cf")

y_limits <- range(df$len, na.rm = TRUE)
pad <- 0.05 * diff(y_limits)
if (!is.finite(pad) || pad == 0) pad <- 0.05 * max(abs(y_limits), 0.01, na.rm = TRUE)
y_limits <- c(y_limits[1] - pad, y_limits[2] + pad)

# Match boxplot dodge and point jitter-dodge (same width + preserve)
dodge_w <- 0.75
pd <- position_dodge(width = dodge_w)

plot_gene <- function(g, show_y_axis = TRUE) {
  ggplot(df %>% filter(gene == g), aes(x = site_type, y = len)) +
    geom_boxplot(
      aes(fill = season),
      outlier.shape = NA,
      position = pd
    ) +
    geom_point(
      aes(shape = site_id, color = season, group = season),
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = dodge_w
      ),
      size = 2,
      stroke = 1
    ) +
    scale_shape_manual(values = c("Dairy" = 5, "Swine & Dairy" = 1, "UI" = 6, "WW" = 4)) +
    scale_fill_manual(values = season_colors, name = "Season") +
    scale_color_manual(values = season_colors, name = "Season") +
    guides(
      color = guide_legend(order = 1),
      fill = guide_none()
    ) +
    stat_compare_means(
      aes(group = season),
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE
    ) +
    coord_cartesian(ylim = y_limits) +
    scale_y_continuous(
      name = if (show_y_axis) "ARG Relative Abundance" else NULL,
      breaks = waiver()
    ) +
    labs(x = NULL, title = g) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = if (show_y_axis) element_text(size = 10) else element_blank(),
      axis.title.y = if (show_y_axis) element_text(size = 12) else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.line.y.left = if (show_y_axis) element_line(color = "black") else element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

genes <- unique(df$gene)
plots <- lapply(seq_along(genes), function(i) {
  plot_gene(genes[i], show_y_axis = (i == 1))
})

wrap_plots(plots, ncol = 3) +
  plot_layout(guides = "collect")
