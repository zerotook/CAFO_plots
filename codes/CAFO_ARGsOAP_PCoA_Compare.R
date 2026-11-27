library(readr)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggrepel)

infile <- "inputs/Compare_combine_merged_with_gene_c.csv"

# 1) Read + clean
arg <- read_csv(infile, show_col_types = FALSE) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

gene_ids <- arg$gene
mat <- arg %>% select(-gene) %>% as.matrix()
rownames(mat) <- gene_ids
mat <- mat[rowSums(mat) > 0, , drop = FALSE]

# 2) Assign site types from sample names (fix D vs sediment_SD... ordering)
sample_names <- colnames(mat)

site_type <- case_when(
  grepl("^sediment",     sample_names, ignore.case = TRUE) ~ "sediment",
  grepl("^surfacewater", sample_names, ignore.case = TRUE) ~ "surfacewater",
  grepl("^SD\\d+$",      sample_names, ignore.case = TRUE) ~ "Swine & Dairy",  # SD1..SD8
  grepl("^D\\d+$",       sample_names, ignore.case = TRUE) ~ "Dairy",          # D1..D8
  grepl("^Reference",      sample_names, ignore.case = TRUE) ~ "Reference",
  grepl("^acs",      sample_names, ignore.case = TRUE) ~ "activated sludge",
  #grepl("^OLW",      sample_names, ignore.case = TRUE) ~ "OLW",
  TRUE ~ NA_character_
)

# Keep only the groups we care about
keep <- !is.na(site_type)
mat <- mat[, keep, drop = FALSE]
site_type <- factor(site_type[keep], levels = c("Dairy", "Swine & Dairy", "surfacewater", "sediment", "Reference", "activated sludge"))#,, "activated sludge", "OLW","hospital swab"
sample_names <- sample_names[keep]

# 3) Bray distance + PCoA
dist_mat <- vegdist(t(mat), method = "bray")
pcoa_res <- cmdscale(dist_mat, k = 2, eig = TRUE)
var_exp <- 100 * pcoa_res$eig / sum(pcoa_res$eig)

pcoa_df <- data.frame(
  Sample   = sample_names,
  PCoA1    = pcoa_res$points[, 1],
  PCoA2    = pcoa_res$points[, 2],
  SiteType = site_type
)

# 4) PERMANOVA label (use SiteType consistently)
set.seed(123)
perm_res <- adonis2(dist_mat ~ SiteType, data = pcoa_df, permutations = 999)
r2 <- perm_res$R2[1]
pval <- perm_res$`Pr(>F)`[1]
label_txt <- paste0("PERMANOVA: R²=", round(r2, 3), ", p=", signif(pval, 3))
print(perm_res)

# 5) Plot + label each sample
# Custom colors
site_colors <- c(
  "Dairy"        = "#a8d3a0",
  "Swine & Dairy"= "#f8b072",
  "surfacewater" = "darkblue",
  "sediment"     = "#a77260",
  "Reference"      = "#b8a88f",
 # "OLW" = "yellow",
  "activated sludge" = "darkred"
)

# Custom shapes (pick any set you like)
site_shapes <- c(
  "Dairy"        = 16,  # filled circle
  "Swine & Dairy"= 17,  # triangle
  "surfacewater" = 15,  # square
  "sediment"     = 18,  # diamond
  "Reference"      = 1,   # star
 # "OLW" = 4,
  "activated sludge"      = 17
)

p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2,
                         color = SiteType, shape = SiteType)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = subset(pcoa_df, SiteType == "Reference"),
    aes(label = Sample),
    size = 3,
    force = 2,
    box.padding = 0.5,
    point.padding = 0.4,
    show.legend = FALSE
  ) +
  scale_color_manual(values = site_colors) +
  scale_shape_manual(values = site_shapes) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_txt,
    hjust = 1.05, vjust = 1.3,
    size = 3.5,
    color = "black"
  ) +
  theme_minimal() +
  labs(
    title = "PCoA of ARG gene profiles (Bray–Curtis)",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    color = "Site type",
    shape = "Site type"
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )

print(p)
