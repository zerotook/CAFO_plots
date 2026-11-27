library(readr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# load your input file here
genera <- read_csv(
  "inputs/Specis_clade_clean.csv",
  name_repair = "unique"
)
# Top 30 overall
abund_top30 <- genera %>%
  select(-lineage) %>%
  mutate(total_abundance = rowSums(across(-name))) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:10) %>%
  select(-total_abundance)

# Ensure unique, explicit row labels
abund_top30 <- abund_top30 %>%
  mutate(name = make.unique(as.character(name)))

row_labels <- abund_top30$name
mat <- abund_top30 %>%
  select(-name) %>%
  as.matrix()

rownames(mat) <- abund_top30$name

# Apply log10(x+1) transform, could adjust to log10(x) by removing +1
#mat_log <- log10(mat+1)

clean_col_labels <- gsub("\\.\\.\\.[0-9]+$", "", colnames(mat))

# Define a color palette
my_colors <- colorRampPalette(c("#5EAFA6", "white", "#DEC27B"))(100)

#adjust cluster rows or columns option by change TRUE or FALSE
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_col = clean_col_labels,
  fontsize_row = 10,
  fontsize_col = 9,
  color = my_colors
)

# If labels still don’t appear in the RStudio Plots pane, save a bigger file:
# pheatmap(
#   mat,
#   scale = "row",
#   cluster_rows = TRUE,
#   cluster_cols = FALSE,
#   show_rownames = TRUE,
#   labels_row = row_labels,
#   fontsize_row = 10,
#   fontsize_col = 9,
#   filename = "Top30_Genera_heatmap.png",  # <-- larger canvas
#   width = 9, height = 12                  # adjust as needed
# )

