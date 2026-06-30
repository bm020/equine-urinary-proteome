# ============================================================
# Figure 2C: Candidate protein heatmap for EGS-increased proteins
# ============================================================
#
# Purpose:
# Generate the Figure 2C heatmap showing relative abundance of
# candidate proteins highlighted from the EGS-increased group.
#
# Input file:
#   - Filtered_final.tsv
#
# Correct biological mapping:
#   U10 = EGS5, excluded
#   U9  = EGS4, excluded
#   U8  = EGS3, retained
#   U7  = EGS2, retained
#   U6  = EGS1, retained
#   U5  = CTRL1
#   U4  = CTRL2
#   U3  = CTRL3
#   U2  = CTRL4
#   U1  = CTRL5
#
# Retained samples:
#   CTRL1, CTRL2, CTRL3, CTRL4, CTRL5, EGS1, EGS2, EGS3
#
# Processing:
#   - Protein quantities are log2-transformed using log2(x + 1).
#   - Sample-wise median normalisation is applied.
#   - Missing PG.Genes values are replaced with Protein_ID.
#   - Candidate proteins are displayed using curated labels.
#   - Final protein order is data-driven, based on decreasing
#     median log2 fold change in EGS versus CTRL urine.
#   - Row z-score scaling is applied for heatmap display.
#   - Z-scores are capped at -2 and 2.
#
# Package/method references:
#   - readr::read_tsv:
#     https://readr.tidyverse.org/reference/read_delim.html
#   - dplyr data manipulation:
#     https://dplyr.tidyverse.org/
#   - tidyr data reshaping:
#     https://tidyr.tidyverse.org/
#   - ggplot2 figure generation:
#     https://ggplot2.tidyverse.org/
#   - svglite SVG export:
#     https://svglite.r-lib.org/
#   - ComplexHeatmap heatmap visualisation:
#     https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
#   - circlize colour mapping:
#     https://cran.r-project.org/package=circlize
#   - row z-score scaling:
#     https://stat.ethz.ch/R-manual/R-devel/library/base/html/scale.html
#
# Output files:
#   - Figure_2C_candidate_heatmap_EGS_increased.png
#   - Figure_2C_candidate_heatmap_EGS_increased.pdf
#   - Figure_2C_candidate_heatmap_EGS_increased.svg
#   - Figure_2C_candidate_heatmap_EGS_increased_selected_proteins.csv
#   - Figure_2C_candidate_heatmap_EGS_increased_mapping.csv
#   - Figure_2C_candidate_heatmap_EGS_increased_sessionInfo.txt
# ============================================================


# ----------------------------
# 0. Package checks
# ----------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "ComplexHeatmap",
  "circlize",
  "grid",
  "svglite"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following packages are missing: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script.\n",
    "For ComplexHeatmap, use: BiocManager::install('ComplexHeatmap')"
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)


# ----------------------------
# 1. User settings
# ----------------------------

input_file <- "Filtered_final.tsv"
output_prefix <- "Figure_2C_candidate_heatmap_EGS_increased"

ctrl_keep <- c("CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5")
egs_keep  <- c("EGS1", "EGS2", "EGS3")
keep_samples <- c(ctrl_keep, egs_keep)

# Candidate protein panel for Figure 2C.
# The final heatmap order is calculated later by decreasing
# median log2FC in EGS versus CTRL urine.
protein_panel <- c(
  "GC",
  "SERPINC1",
  "Pepsin A",
  "ALB",
  "HP",
  "AHSG",
  "TF",
  "CTRC",
  "SERPINA3",
  "HBB",
  "MSMB",
  "NPC2",
  "VMO1",
  "LMAN2",
  "Q95182",
  "LCN2",
  "DSC3",
  "CFB",
  "LRG1",
  "AZGP1"
)

group_cols <- c(
  "CTRL" = "#15a2a2",
  "EGS"  = "#fea802"
)

heat_cols <- colorRamp2(
  c(-2, 0, 2),
  c("#0072B2", "white", "#D55E00")
)


# ----------------------------
# 2. Read filtered protein table
# ----------------------------

stopifnot(file.exists(input_file))

Filtered_final <- read_tsv(input_file, show_col_types = FALSE)


# ----------------------------
# 3. Detect PG.Quantity columns
# ----------------------------

sample_cols <- grep(
  "Horse_urine_U.*PG\\.Quantity$",
  colnames(Filtered_final),
  value = TRUE
)

stopifnot(length(sample_cols) == 10)

cat("\nDetected PG.Quantity columns:\n")
print(sample_cols)


# ----------------------------
# 4. Rename columns using correct biological mapping
# ----------------------------

mapping_tbl <- tibble(
  original_col = sample_cols,
  new_name = c(
    "EGS5_excluded",  # U10
    "EGS4_excluded",  # U9
    "EGS3",           # U8
    "EGS2",           # U7
    "EGS1",           # U6
    "CTRL1",          # U5
    "CTRL2",          # U4
    "CTRL3",          # U3
    "CTRL4",          # U2
    "CTRL5"           # U1
  )
)

cat("\nSample mapping used for Figure 2C:\n")
print(mapping_tbl)

colnames(Filtered_final)[match(sample_cols, colnames(Filtered_final))] <- mapping_tbl$new_name

stopifnot(all(keep_samples %in% colnames(Filtered_final)))


# ----------------------------
# 5. Build wide protein table
# ----------------------------

data_wide <- Filtered_final %>%
  select(
    Protein_ID = PG.ProteinGroups,
    Gene = PG.Genes,
    Protein_description = PG.ProteinDescriptions,
    all_of(keep_samples)
  )


# ----------------------------
# 6. Clean gene symbols
# ----------------------------
# Q95182 has missing PG.Genes in the dataset, so missing Gene
# values are replaced with Protein_ID to prevent loss of this entry.
# ----------------------------

data_wide <- data_wide %>%
  mutate(
    Protein_ID = as.character(Protein_ID),
    Gene = as.character(Gene),
    Protein_description = as.character(Protein_description),
    Gene = str_replace_all(Gene, ";", ","),
    Gene = if_else(
      is.na(Gene) | Gene == "",
      Protein_ID,
      Gene
    )
  ) %>%
  separate_rows(Gene, sep = ",") %>%
  mutate(Gene = str_trim(Gene)) %>%
  filter(!is.na(Gene), Gene != "")


# ----------------------------
# 7. Convert to long format and log2-transform
# ----------------------------

data_long <- data_wide %>%
  pivot_longer(
    cols = all_of(keep_samples),
    names_to = "Sample",
    values_to = "Quantity"
  ) %>%
  mutate(
    Quantity = as.numeric(Quantity),
    log2_abundance = log2(Quantity + 1),
    Group = case_when(
      Sample %in% ctrl_keep ~ "CTRL",
      Sample %in% egs_keep  ~ "EGS",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))


# ----------------------------
# 8. Median-normalise within each sample
# ----------------------------

data_long <- data_long %>%
  group_by(Sample) %>%
  mutate(
    median_normalised_log2 =
      log2_abundance - median(log2_abundance, na.rm = TRUE)
  ) %>%
  ungroup()


# ----------------------------
# 9. Collapse duplicate protein/gene entries per sample
# ----------------------------

safe_max <- function(x) {
  if (all(is.na(x) | !is.finite(x))) {
    NA_real_
  } else {
    max(x, na.rm = TRUE)
  }
}

data_gene <- data_long %>%
  group_by(Protein_ID, Gene, Protein_description, Sample, Group) %>%
  summarise(
    abundance = safe_max(median_normalised_log2),
    .groups = "drop"
  ) %>%
  filter(is.finite(abundance))


# ----------------------------
# 10. Create curated display labels
# ----------------------------

data_gene <- data_gene %>%
  mutate(
    Gene_display = case_when(
      Gene == "LOC100062179" ~ "Pepsin A",
      Gene == "LOC100067869" ~ "HP",
      Protein_ID == "Q95182" ~ "Q95182",
      Gene == "Q95182" ~ "Q95182",
      str_detect(
        Protein_description,
        regex("Major allergen Equ c 1|Equ c 1", ignore_case = TRUE)
      ) ~ "Q95182",
      TRUE ~ Gene
    )
  )

data_gene_display <- data_gene %>%
  group_by(Gene_display, Sample, Group) %>%
  summarise(
    abundance = safe_max(abundance),
    .groups = "drop"
  ) %>%
  filter(is.finite(abundance))


# ----------------------------
# 11. Check selected protein panel
# ----------------------------

available_labels <- unique(data_gene_display$Gene_display)

missing_panel <- setdiff(protein_panel, available_labels)

if (length(missing_panel) > 0) {
  warning(
    "The following requested protein labels were not found in the filtered data: ",
    paste(missing_panel, collapse = ", ")
  )
}

protein_panel_found <- protein_panel[protein_panel %in% available_labels]

if (length(protein_panel_found) == 0) {
  stop("None of the requested protein panel labels were found in the data.")
}

cat("\nRequested protein panel found in data:\n")
print(protein_panel_found)


# ----------------------------
# 12. Calculate data-driven protein order
# ----------------------------

selected_protein_table <- data_gene_display %>%
  filter(Gene_display %in% protein_panel_found) %>%
  group_by(Gene_display, Group) %>%
  summarise(
    median_abundance = median(abundance, na.rm = TRUE),
    n_detected = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(median_abundance, n_detected)
  ) %>%
  mutate(
    log2FC_EGS_vs_CTRL = median_abundance_EGS - median_abundance_CTRL,
    FC_EGS_vs_CTRL = 2^log2FC_EGS_vs_CTRL
  ) %>%
  arrange(desc(log2FC_EGS_vs_CTRL))

cat("\nSelected protein table ordered by decreasing median log2FC:\n")
print(selected_protein_table)

write_csv(
  selected_protein_table,
  paste0(output_prefix, "_selected_proteins.csv")
)

protein_panel_ordered <- selected_protein_table %>%
  filter(!is.na(log2FC_EGS_vs_CTRL)) %>%
  arrange(desc(log2FC_EGS_vs_CTRL)) %>%
  pull(Gene_display)

cat("\nFinal heatmap protein order:\n")
print(protein_panel_ordered)


# ----------------------------
# 13. Save accession/gene/display-label mapping
# ----------------------------

selected_mapping_table <- data_gene %>%
  filter(Gene_display %in% protein_panel_ordered) %>%
  select(
    Protein_ID,
    Gene,
    Gene_display,
    Protein_description
  ) %>%
  distinct() %>%
  arrange(match(Gene_display, protein_panel_ordered))

cat("\nSelected protein mapping table:\n")
print(selected_mapping_table)

write_csv(
  selected_mapping_table,
  paste0(output_prefix, "_mapping.csv")
)


# ----------------------------
# 14. Create heatmap matrix
# ----------------------------

heatmap_mat <- data_gene_display %>%
  filter(Gene_display %in% protein_panel_ordered) %>%
  select(Gene_display, Sample, abundance) %>%
  pivot_wider(
    names_from = Sample,
    values_from = abundance
  ) %>%
  as.data.frame()

rownames(heatmap_mat) <- heatmap_mat$Gene_display
heatmap_mat$Gene_display <- NULL

heatmap_mat <- as.matrix(heatmap_mat)

# Force CTRL then EGS sample order
heatmap_mat <- heatmap_mat[, keep_samples, drop = FALSE]

# Force data-driven row order
heatmap_mat <- heatmap_mat[protein_panel_ordered, , drop = FALSE]


# ----------------------------
# 15. Row z-score scaling
# ----------------------------

row_scale <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  s <- sd(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / s
}

heatmap_scaled <- t(apply(heatmap_mat, 1, row_scale))

# Cap extreme values for cleaner visualisation
heatmap_scaled[heatmap_scaled > 2]  <- 2
heatmap_scaled[heatmap_scaled < -2] <- -2

# Replace remaining NA with 0 for display only
heatmap_scaled[is.na(heatmap_scaled)] <- 0


# ----------------------------
# 16. Sample annotation
# ----------------------------

sample_group <- factor(
  c(rep("CTRL", length(ctrl_keep)), rep("EGS", length(egs_keep))),
  levels = c("CTRL", "EGS")
)

names(sample_group) <- keep_samples

top_ha <- HeatmapAnnotation(
  Group = sample_group,
  col = list(Group = group_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = unit(4, "mm"),
  annotation_height = unit(4, "mm")
)


# ----------------------------
# 17. Legends
# ----------------------------

lgd_relative <- Legend(
  title = "Relative\nabundance",
  col_fun = heat_cols,
  at = c(-2, -1, 0, 1, 2),
  labels = c("-2", "-1", "0", "1", "2"),
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  legend_height = unit(34, "mm"),
  grid_width = unit(4, "mm")
)

lgd_group <- Legend(
  title = "Group",
  labels = c("CTRL", "EGS"),
  legend_gp = gpar(fill = group_cols[c("CTRL", "EGS")]),
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  grid_height = unit(4, "mm"),
  grid_width = unit(4, "mm")
)


# ----------------------------
# 18. Heatmap object
# ----------------------------

ht <- Heatmap(
  heatmap_scaled,
  name = "Relative abundance",
  col = heat_cols,
  top_annotation = top_ha,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  column_split = sample_group,
  column_gap = unit(4, "mm"),

  column_title = " ",
  column_title_gp = gpar(fontsize = 0),

  show_row_dend = FALSE,
  show_column_dend = FALSE,

  row_names_side = "right",
  row_names_gp = gpar(
    fontsize = 7.8,
    fontface = "plain",
    col = "black"
  ),
  row_names_max_width = unit(42, "mm"),

  column_names_gp = gpar(
    fontsize = 8.5,
    fontface = "plain",
    col = "black"
  ),
  column_names_rot = 45,
  column_names_centered = TRUE,

  show_heatmap_legend = FALSE,

  border = TRUE,
  rect_gp = gpar(
    col = "grey88",
    lwd = 0.35
  ),

  row_title = NULL,

  width = unit(92, "mm"),
  height = unit(112, "mm")
)


# ----------------------------
# 19. Draw function
# ----------------------------

draw_candidate_heatmap <- function() {
  draw(
    ht,
    heatmap_legend_list = list(lgd_relative, lgd_group),
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    align_heatmap_legend = "heatmap_top",
    align_annotation_legend = "heatmap_top",
    merge_legends = TRUE,
    legend_grouping = "original",
    padding = unit(c(1, 1, 1, 1), "mm")
  )
}


# ----------------------------
# 20. Preview
# ----------------------------

grid.newpage()
draw_candidate_heatmap()


# ----------------------------
# 21. Export PNG, PDF and SVG
# ----------------------------

png(
  filename = paste0(output_prefix, ".png"),
  width = 7.1,
  height = 6.0,
  units = "in",
  res = 600
)

grid.newpage()
draw_candidate_heatmap()
dev.off()

pdf(
  file = paste0(output_prefix, ".pdf"),
  width = 7.1,
  height = 6.0,
  useDingbats = FALSE
)

grid.newpage()
draw_candidate_heatmap()
dev.off()

svglite::svglite(
  filename = paste0(output_prefix, ".svg"),
  width = 7.1,
  height = 6.0
)

grid.newpage()
draw_candidate_heatmap()
dev.off()


# ----------------------------
# 22. Save session information
# ----------------------------

sink(paste0(output_prefix, "_sessionInfo.txt"))
sessionInfo()
sink()


cat("\nSaved:\n")
cat("-", paste0(output_prefix, ".png"), "\n")
cat("-", paste0(output_prefix, ".pdf"), "\n")
cat("-", paste0(output_prefix, ".svg"), "\n")
cat("-", paste0(output_prefix, "_selected_proteins.csv"), "\n")
cat("-", paste0(output_prefix, "_mapping.csv"), "\n")
cat("-", paste0(output_prefix, "_sessionInfo.txt"), "\n")
