library(VennDiagram)
library(grid)

df <- read.table(
  "inputs/combined_normalized_cell.subtype.tsv",
  header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)

# Force numeric for all sample columns (everything except subtype)
df[-1] <- lapply(df[-1], function(x) as.numeric(x))
df[is.na(df)] <- 0

# Define sample groups by ID (number before F/S)
dairy_ids       <- c("1", "2", "3", "4", "7")
swine_dairy_ids <- c("5", "6", "8", "9")
clean_ids       <- c("10")
waste_ids       <- c("11")

get_cols_for_ids <- function(ids, all_names) {
  # Extract leading digits as site id: "10S" -> "10", "1F" -> "1"
  prefix <- sub("^([0-9]+).*", "\\1", all_names)
  all_names[prefix %in% ids]
}

all_cols <- colnames(df)[-1]

dairy_cols       <- get_cols_for_ids(dairy_ids, all_cols)
swine_dairy_cols <- get_cols_for_ids(swine_dairy_ids, all_cols)
clean_cols       <- get_cols_for_ids(clean_ids, all_cols)
waste_cols       <- get_cols_for_ids(waste_ids, all_cols)

# Sanity check
print(list(
  dairy = dairy_cols,
  swine_dairy = swine_dairy_cols,
  clean = clean_cols,
  waste = waste_cols
))

present_in_group <- function(colset) {
  if (length(colset) == 0) return(character(0))
  sub_df <- df[, colset, drop = FALSE]
  df$subtype[rowSums(sub_df > 0) > 0]
}

set_dairy       <- present_in_group(dairy_cols)
set_swine_dairy <- present_in_group(swine_dairy_cols)
set_clean       <- present_in_group(clean_cols)
set_waste       <- present_in_group(waste_cols)

# ---- Save to PDF (recommended) ----
pdf("examples/CAFO_Venn_v2_example.pdf", width = 8, height = 8)

venn.plot <- venn.diagram(
  x = list(
    Dairy             = set_dairy,
    "Swine & Dairy"   = set_swine_dairy,
    UI                = set_clean,
    "ww"              = set_waste
  ),
  filename = NULL,
  fill = c("#a8d3a0", "#f8b072", "#4f99c9", "#ec3e31"),
  alpha = 0.8,
  cex = 1.5,
  fontfamily = "sans",   # intersection font family
  cat.cex = 1,          # category label font size
  cat.fontfamily = "sans", # category label font
  cat.fontface = "bold",
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)
dev.off()
