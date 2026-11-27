suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(circlize)
})

circos.clear()

# ===================== USER SETTINGS =====================
infile    <- "inputs/combined_D_spring_cleaned.tsv"
tax_col   <- "phylum"   # use "phylum_ncbi" for phylum; change if needed

# =========================================================

# ---------- FIXED COLOR MAPS (stable across files) ----------
# ARG classes (multidrug subtypes will be collapsed to "multidrug")
ARG_COLORS <- c(
  "beta-lactam"     = "#FF9896",
  "macrolide-lincosamide-streptogramin"             = "#9467BD",
  "defensin"        = "#2CA02C",
  "polymyxin"       = "#FF7F0E",
  "rifamycin"       = "#8C564B",
  "fosfomycin"      = "#E377C2",
  "tetracycline"    = "#D62728",
  "trimethoprim"    = "#BCBD22",
  "vancomycin"      = "#17BECF",
  "aminoglycoside"  = "#AEC7E8",
  "chloramphenicol" = "#98DF8A",
  "sulfonamide"     = "#F5A834",
  "quinolone"       = "#C5B0D5",
  "multidrug"       = "#EBE9E8",
  "bacitracin"      = "#8DD3C7",
  "mupirocin"       = "#F7B6D2", 
  "novobiocin"      = "#DBDB8D"
)

# Phylum colors
# Phylum colors (traditional names)
PHYLUM_COLORS <- c(
  "Proteobacteria"        = "#4E79A7",
  "Bacteroidetes"         = "#F28E2B",
  "Bdellovibrio"          = "#E15759",
  "Firmicutes"            = "#76B7B2",
  "Planctomycetes"        = "#59A14F",
  "Thermodesulfobacteria" = "#EDC948",
  "Myxobacteria"          = "#B07AA1",
  "Cyanobacteria"         = "purple",
  "Epsilonproteobacteria" = "#9C755F",
  "Actinobacteria"        = "darkgreen",
  "Acidobacteria"         = "#BAB0AC",
  "Fusobacteria"          = "#86BCB6",
  "Chordata"              = "#D4A6C8",
  "Verrucomicrobia"       = "#8CD17D",
  "Chlamydiae"            = "#B6992D",
  "Chlorobi"              = "#499894",
  "Unclassified"          = "#BDBDBD"
)


FALLBACK_COLOR <- "#D9D9D9"  # used if a new/unknown category appears
# ------------------------------------------------------------

# Display name mapping for ARG classes (for plot labels only)
ARG_DISPLAY_NAMES <- c(
  "macrolide-lincosamide-streptogramin" = "MLS"
  # Add more mappings here if needed in the future
  # "another-long-name" = "Short"
)

# Helper function to get display name (returns original if no mapping exists)
get_display_name <- function(name) {
  if (name %in% names(ARG_DISPLAY_NAMES)) {
    return(ARG_DISPLAY_NAMES[name])
  }
  return(name)
}

# Helper to find a column by name, case/underscore-insensitive
find_col <- function(nms, target){
  idx <- which(tolower(gsub("[ _]", "", nms)) == gsub("[ _]", "", tolower(target)))
  if(length(idx) == 1) return(nms[idx])
  variants <- tolower(c(target, paste0("tax_", target)))
  idx2 <- which(tolower(gsub("[ _]", "", nms)) %in% gsub("[ _]", "", variants))
  if(length(idx2) >= 1) return(nms[idx2[1]])
  NA_character_
}

# Make sure every level has a color (stable: known ones fixed, unknowns fallback)
ensure_colors <- function(levels_vec, fixed_map, fallback = FALLBACK_COLOR){
  out <- rep(fallback, length(levels_vec))
  names(out) <- levels_vec
  hit <- intersect(levels_vec, names(fixed_map))
  out[hit] <- fixed_map[hit]
  out
}

# ---- Read ----
edges <- read_tsv(infile, show_col_types = FALSE) %>%
  rename(arg_class = to, tpm = value)

taxon_colname <- find_col(names(edges), tax_col)
stopifnot(!is.na(taxon_colname))

edges <- edges %>%
  mutate(
    arg_class = trimws(as.character(arg_class)),
    tpm       = suppressWarnings(as.numeric(tpm))
  ) %>%
  rename(taxon = !!taxon_colname) %>%
  mutate(taxon = trimws(as.character(taxon)))

# ---- Map GTDB / -ota phylum names to traditional names ----
phylum_map <- c(
  "Pseudomonadota"          = "Proteobacteria",
  "Bacteroidota"            = "Bacteroidetes",
  "Bdellovibrionota"        = "Bdellovibrio",
  "Bacillota"               = "Firmicutes",
  "Planctomycetota"         = "Planctomycetes",
  "Thermodesulfobacteriota" = "Thermodesulfobacteria",
  "Myxococcota"             = "Myxobacteria",
  "Cyanobacteriota"         = "Cyanobacteria",
  "Campylobacterota"        = "Epsilonproteobacteria",
  "Actinomycetota"          = "Actinobacteria",
  "Acidobacteriota"         = "Acidobacteria",
  "Fusobacteriota"          = "Fusobacteria",
  "Chordata"                = "Chordata",
  "Verrucomicrobiota"       = "Verrucomicrobia",
  "Chlamydiota"             = "Chlamydiae",
  "Chlorobiota"             = "Chlorobi",
  "unclassified"            = "Unclassified"
)

edges <- edges %>%
  mutate(
    taxon = ifelse(taxon %in% names(phylum_map),
                   phylum_map[taxon],
                   taxon)
  )


stopifnot(all(c("arg_class","taxon","tpm") %in% names(edges)))

# ---- Combine multidrug subclasses into one ----
# e.g., multidrug@ABC, multidrug@MATE, multidrug@MFS, multidrug@SMR -> multidrug
edges <- edges %>%
  mutate(arg_class = ifelse(grepl("^multidrug@", arg_class), "multidrug", arg_class))

# ---- Filter low TPM BEFORE aggregation ----
tpm_cut <- quantile(edges$tpm, 0, na.rm = TRUE)
edges_f <- edges %>% filter(!is.na(tpm), tpm > tpm_cut)

# ---- Aggregate ----
edges_sum <- edges_f %>%
  group_by(arg_class, taxon) %>%
  summarise(weight = sum(tpm, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    arg_class = trimws(as.character(arg_class)),
    taxon     = trimws(as.character(taxon))
  ) %>%
  filter(!is.na(arg_class), arg_class != "", !is.na(taxon), taxon != "")

# ---- Normalize (0–1) ----
max_w <- max(edges_sum$weight, na.rm = TRUE)
edges_sum <- edges_sum %>%
  mutate(TPM_scaled = ifelse(is.finite(max_w) && max_w > 0, weight / max_w, 0))

# ---- Sector order ----
args <- sort(unique(edges_sum$arg_class))
taxa <- sort(unique(edges_sum$taxon))
order <- c(args, taxa)

all_sectors <- sort(unique(c(edges_sum$arg_class, edges_sum$taxon)))
missing <- setdiff(all_sectors, order)
if (length(missing) > 0) stop(sprintf("Sectors missing from `order`: %s", paste(missing, collapse = ", ")))

# ---- Colors applied (stable across files) ----
arg_cols   <- ensure_colors(args, ARG_COLORS)
tax_cols   <- ensure_colors(taxa, PHYLUM_COLORS)
grid_cols  <- c(arg_cols, tax_cols)

# Color links by ARG class (recommended)
link_cols <- unname(arg_cols[edges_sum$arg_class])

# ---- Gaps ----
gap_after <- c(
  if (length(args) > 0) rep(2, max(0, length(args) - 1)) else numeric(0),
  10,
  if (length(taxa) > 0) rep(2, max(0, length(taxa) - 1)) else numeric(0),
  10
)

# ---- Plot ----
circos.clear()
circos.par(
  start.degree = 120,
  gap.after    = gap_after,
  track.margin = c(0.02, 0.02)
)

chordDiagram(
  x = edges_sum[, c("arg_class", "taxon", "TPM_scaled")],
  order = order,
  grid.col = grid_cols,     # fixed sector colors
  col = link_cols,          # links colored by ARG class
  transparency = 0.25,
  annotationTrack = "grid",
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  link.arr.type = "big.arrow",
  preAllocateTracks = list(track.height = 0.10)
)

# ---- Labels ----
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  track.height = 0.15,
  panel.fun = function(x, y) {
    sector.name <- CELL_META$sector.index
    # Use display name for ARG classes, original name for taxa
    display.name <- if (sector.name %in% args) {
      get_display_name(sector.name)
    } else {
      sector.name
    }
    
    if (sector.name %in% args) {
      circos.text(
        CELL_META$xcenter, CELL_META$ycenter, display.name,
        facing = "downward", niceFacing = TRUE, adj = c(0.1, 0.2), cex = 0.75
      )
    } else {
      xlim <- CELL_META$xlim; ylim <- CELL_META$ylim
      circos.text(
        mean(xlim), mean(ylim), display.name,
        facing = "downward", niceFacing = TRUE, adj = c(0.85, 0.8), cex = 0.75
      )
    }
  }
)

title("Dairy Spring")
#750*600


# ---- Optional: print the color keys used in THIS plot ----
#message("\nARG colors used:")
#print(arg_cols)
#message("\nTaxon colors used:")
#print(tax_cols)
