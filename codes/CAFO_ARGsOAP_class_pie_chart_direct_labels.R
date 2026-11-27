# install.packages(c("circlize", "readr", "dplyr", "RColorBrewer"))
library(circlize)
library(readr)
library(dplyr)
library(RColorBrewer)

input_file <- "inputs/combined_normalized_cell.subtype_circle - 10_11.csv"
output_dir <- "examples"

output_png <- file.path(output_dir, "ARG_class_pie_chart_direct_labels.png")
output_all_classes_csv <- file.path(output_dir, "ARG_class_totals_all_classes.csv")
output_pie_csv <- file.path(output_dir, "ARG_class_pie_data_top15.csv")

# Same class-level filter as the ring plot: after removing excluded classes,
# keep the top 15 ARG classes by summed abundance.
top_classes <- 15

exclude_classes <- c("biocide", "defensin", "ionophore")

# Same text-style control used in your original circos.text() code.
type_label_cex <- 1.8

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

cat("Reading input file...\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(input_file, show_col_types = FALSE) %>%
  rename(subtype_full = 1)

sample_cols <- setdiff(names(df), "subtype_full")

dat <- df %>%
  mutate(
    arg_class = sub("__.*$", "", subtype_full),
    value = rowSums(across(all_of(sample_cols)), na.rm = TRUE)
  ) %>%
  filter(
    is.finite(value),
    value > 0,
    !arg_class %in% exclude_classes
  )

class_totals <- dat %>%
  group_by(arg_class) %>%
  summarise(total = sum(value), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(percent = total / sum(total))

write_csv(class_totals, output_all_classes_csv)

pie_data <- class_totals %>%
  head(min(top_classes, nrow(class_totals))) %>%
  mutate(percent = total / sum(total))

write_csv(pie_data, output_pie_csv)

classes <- pie_data$arg_class
xlim <- cbind(rep(0, length(classes)), pie_data$total)
class_col <- custom_colors[classes]

missing_idx <- is.na(class_col)
if (any(missing_idx)) {
  missing_classes <- classes[missing_idx]
  cat("Warning: No custom color found for:", paste(missing_classes, collapse = ", "), "\n")

  fallback_pal <- colorRampPalette(brewer.pal(8, "Set2"))(sum(missing_idx))
  class_col[missing_idx] <- fallback_pal
}

names(class_col) <- classes

cat("Generating ARG class pie chart with circos.text labels...\n")

png(output_png, width = 3500, height = 3500, res = 300, units = "px")

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 1.5,
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0, 0),
  points.overflow.warning = FALSE
)

circos.initialize(factors = classes, xlim = xlim)

circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.98,
  bg.border = NA,
  panel.fun = function(x, y) {
    t <- CELL_META$sector.index
    xlim <- CELL_META$xlim

    circos.rect(
      xlim[1], 0,
      xlim[2], 1,
      col = class_col[[t]],
      border = "white",
      lwd = 1.2
    )

    circos.text(
      mean(xlim), 0.77,
      t,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0.5, 0.5),
      cex = type_label_cex,
      col = "black",
      font = 2
    )
  }
)

mtext(
  "ARG Class Distribution",
  side = 3,
  line = -1,
  cex = 1.4,
  font = 2
)

circos.clear()
dev.off()

cat("Saved pie chart as:", output_png, "\n")
cat("Saved all class totals as:", output_all_classes_csv, "\n")
cat("Saved plotted pie data as:", output_pie_csv, "\n")
