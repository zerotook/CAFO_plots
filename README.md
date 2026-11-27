# CAFO plotting bundle

This folder brings together the 12 R plotting scripts from `final codes` and the 11 input files they reference. The original scripts and inputs remain in their original locations.

## Folder layout

- `codes/` — copies of the R scripts. Their input paths now point to `inputs/` relative to the bundle root.
- `inputs/` — copies of the CSV, TSV, and Excel source files.
- `examples/` — example figure outputs (PDF and PNG) generated from the bundled scripts and inputs, plus the pie-chart summary tables. Multi-page PDFs include each panel/page from the original output.

Set the bundle root as your R working directory before running a script so its relative input path resolves. The `qplot.xlsx` script reads the first worksheet.

## Script to input map

| Script | Input file |
| --- | --- |
| `CAFO_ACC_plot_v5.R` | `combined_D_spring_cleaned.tsv` |
| `CAFO_ARGsOAP_class_pie_chart_direct_labels.R` | `combined_normalized_cell.subtype_circle - 10_11.csv` |
| `CAFO_ARGsOAP_PCoA_Compare.R` | `Compare_combine_merged_with_gene_c.csv` |
| `CAFO_ARGsOAP_subtype_barplot.R` | `combined_normalized_cell.subtype_circle - 10_11.csv` |
| `CAFO_ARG_rank1_bubble.R` | `rank_ARG_bubble.csv` |
| `CAFO_ARG_rank1_scatter.R` | `Rank1_scatterplot.csv` |
| `CAFO_ARG_RA_fold_season.R` | `ARGsOAP_RA.csv` |
| `CAFO_BacterialCommunity_heatplot_top30.R` | `Specis_clade_clean.csv` |
| `CAFO_BacterialCommunity_stackedbarchart_clade_v1.R` | `Class_clade_clean.csv` |
| `CAFO_correlation_triangle.R` | `Cross_method.csv` |
| `CAFO_qPCR_boxwhisker_v2.R` | `qplot.xlsx` |
| `CAFO_Venn_v2.R` | `combined_normalized_cell.subtype.tsv` |

Run the scripts with this folder as the working directory. Figure and summary-table outputs are directed to `examples/` where the script has an explicit save step. The example figures are derived from the bundled data; the `inputs/` files themselves remain separate from the examples.

