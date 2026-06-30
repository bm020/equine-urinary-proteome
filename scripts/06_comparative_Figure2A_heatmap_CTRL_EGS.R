# ============================================================
# Figure 2A: Heatmap of CTRL vs EGS urinary proteomic landscape
# ============================================================
#
# Purpose:
# Generate the Figure 2A heatmap showing relative protein
# abundance across five CTRL and three retained EGS urine samples.
#
# Input files:
#   - Filtered_final.tsv
#   - Urine Phase 2 Metadata table.xlsx
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
# Processing:
#   - Retain CTRL1-CTRL5 and EGS1-EGS3 only.
#   - Keep proteins detected in at least 50% of retained samples.
#   - Log2-transform protein quantities using log2(x + 1).
#   - Median-normalise within each sample.
#   - Impute missing values using row medians.
#   - Apply row z-score scaling for heatmap display.
#   - Clamp z-scores to the range -2 to 2.
#
# Output files:
#   - Figure_2A_heatmap_CTRL_EGS.png
#   - Figure_2A_heatmap_CTRL_EGS.pdf
#   - Figure_2A_heatmap_CTRL_EGS.svg
#   - Figure_2A_heatmap_CTRL_EGS_input_matrix.csv
#   - Figure_2A_heatmap_CTRL_EGS_sessionInfo.txt
# ============================================================


# -----------------------------
# 0. Load packages
# -----------------------------

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(pheatmap)
library(viridis)
library(svglite)
library(grid)
library(gtable)


# -----------------------------
# 1. File paths
# -----------------------------

pg_file   <- "Filtered_final.tsv"
meta_file <- "Urine Phase 2 Metadata table.xlsx"

stopifnot(file.exists(pg_file))
stopifnot(file.exists(meta_file))


# -----------------------------
# 2. Read metadata
# -----------------------------

meta_raw <- read_excel(meta_file)

names(meta_raw) <- names(meta_raw) |>
  str_replace_all("\\s+", "_") |>
  str_replace_all("[^A-Za-z0-9_]", "") |>
  str_to_lower()

meta <- meta_raw %>%
  rename(
    sample = any_of(c("sample_id", "sample", "sampleid")),
    age    = any_of(c("age_years", "age", "ageyears")),
    group  = any_of(c("group", "status"))
  ) %>%
  mutate(
    sample = as.character(sample),
    age = suppressWarnings(as.numeric(age)),
    group = case_when(
      str_to_lower(group) %in% c("control", "ctrl") ~ "CTRL",
      str_to_lower(group) %in% c("egs") ~ "EGS",
      TRUE ~ as.character(group)
    )
  )


# -----------------------------
# 3. Read filtered proteomics file
# -----------------------------

pg <- read_tsv(pg_file, show_col_types = FALSE)


# -----------------------------
# 4. Detect and rename PG.Quantity columns
# -----------------------------

sample_cols <- grep(
  "Horse_urine_U.*PG\\.Quantity$",
  colnames(pg),
  value = TRUE
)

stopifnot(length(sample_cols) == 10)

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

cat("\nSample mapping used for Figure 2A:\n")
print(mapping_tbl)

colnames(pg)[match(sample_cols, colnames(pg))] <- mapping_tbl$new_name


# -----------------------------
# 5. Define retained sample order
# -----------------------------

ctrl_keep <- c("CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5")
egs_keep  <- c("EGS1", "EGS2", "EGS3")
keep_samples <- c(ctrl_keep, egs_keep)

stopifnot(all(keep_samples %in% colnames(pg)))

cat("\nRetained samples used in Figure 2A:\n")
print(keep_samples)

cat("\nProtein counts in filtered dataset:\n")
print(colSums(!is.na(pg[, keep_samples]) & is.finite(as.matrix(pg[, keep_samples])) & as.matrix(pg[, keep_samples]) > 0))


# -----------------------------
# 6. Prepare metadata for retained samples
# -----------------------------

meta_use <- meta %>%
  filter(sample %in% keep_samples) %>%
  mutate(
    Group = factor(group, levels = c("CTRL", "EGS")),
    Age = case_when(
      age >= 2  & age <= 5  ~ "2-5",
      age >= 12 & age <= 15 ~ "12-15",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Age)) %>%
  mutate(
    Age = factor(Age, levels = c("2-5", "12-15"))
  )

stopifnot(all(keep_samples %in% meta_use$sample))

meta_use <- meta_use %>%
  mutate(sample = factor(sample, levels = keep_samples)) %>%
  arrange(sample)


# -----------------------------
# 7. Build protein matrix
# -----------------------------

safe_max <- function(x) {
  if (all(is.na(x) | !is.finite(x))) {
    NA_real_
  } else {
    max(x, na.rm = TRUE)
  }
}

mat_tbl <- pg %>%
  select(
    Protein_ID = PG.ProteinGroups,
    all_of(keep_samples)
  ) %>%
  mutate(
    across(all_of(keep_samples), ~ as.numeric(.))
  ) %>%
  group_by(Protein_ID) %>%
  summarise(
    across(all_of(keep_samples), safe_max),
    .groups = "drop"
  )

mat <- mat_tbl %>%
  column_to_rownames("Protein_ID") %>%
  as.matrix()

mode(mat) <- "numeric"


# -----------------------------
# 8. Filter proteins detected in at least 50% of retained samples
# -----------------------------

keep_proteins <- rowSums(!is.na(mat) & is.finite(mat) & mat > 0) >= ceiling(ncol(mat) / 2)

mat <- mat[keep_proteins, , drop = FALSE]

cat("\nProteins retained for Figure 2A heatmap:", nrow(mat), "\n")


# -----------------------------
# 9. Log2 transform and median-normalise
# -----------------------------

mat_log2 <- log2(mat + 1)

sample_medians <- apply(
  mat_log2,
  2,
  function(x) median(x[is.finite(x)], na.rm = TRUE)
)

mat_norm <- sweep(
  mat_log2,
  2,
  sample_medians,
  FUN = "-"
)


# -----------------------------
# 10. Row-median imputation
# -----------------------------

mat_imp <- t(apply(mat_norm, 1, function(v) {
  v[!is.finite(v) | is.na(v)] <- median(v[is.finite(v)], na.rm = TRUE)
  v
}))

colnames(mat_imp) <- colnames(mat_norm)
rownames(mat_imp) <- rownames(mat_norm)


# -----------------------------
# 11. Row z-score transformation for heatmap display
# -----------------------------

mat_z <- t(scale(t(mat_imp)))
mat_z[!is.finite(mat_z)] <- 0

mat_z <- pmax(pmin(mat_z, 2), -2)


# -----------------------------
# 12. Column annotation
# -----------------------------

ann_col <- meta_use %>%
  transmute(
    Sample = as.character(sample),
    Group = Group,
    Age = Age
  ) %>%
  column_to_rownames("Sample")

ann_col <- ann_col[colnames(mat_z), , drop = FALSE]

ann_colors <- list(
  Group = c(
    "CTRL" = "#15a2a2",
    "EGS"  = "#fea802"
  ),
  Age = c(
    "2-5"   = "#D9D9D9",
    "12-15" = "#737373"
  )
)


# -----------------------------
# 13. Add spacing between annotation and heatmap legend
# -----------------------------

add_legend_gap <- function(phm, gap_mm = 8) {
  gt <- phm$gtable

  ann_idx <- which(gt$layout$name %in% c("annotation_legend", "annotation_legend_right"))
  leg_idx <- which(gt$layout$name %in% c("legend", "legend_right"))

  if (length(ann_idx) == 0 || length(leg_idx) == 0) return(phm)

  insert_after_row <- max(gt$layout$b[ann_idx])

  gt2 <- gtable::gtable_add_rows(
    gt,
    heights = grid::unit(gap_mm, "mm"),
    pos = insert_after_row
  )

  phm$gtable <- gt2
  phm
}


# -----------------------------
# 14. Plot heatmap
# -----------------------------

breaks <- seq(-2, 2, length.out = 101)
cols <- viridis(length(breaks) - 1, option = "D")

p_heatmap <- pheatmap(
  mat_z,
  color = cols,
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = length(ctrl_keep),
  annotation_col = ann_col[, c("Age", "Group"), drop = FALSE],
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  border_color = NA,
  fontsize = 9,
  fontsize_col = 9,
  angle_col = 90,
  main = "",
  legend = TRUE,
  treeheight_row = 50,
  treeheight_col = 0,
  silent = TRUE
)

p_heatmap <- add_legend_gap(p_heatmap, gap_mm = 8)


# -----------------------------
# 15. Save helper
# -----------------------------

draw_pheatmap <- function(phm) {
  if (is.list(phm) && !is.null(phm$gtable)) {
    grid::grid.draw(phm$gtable)
  } else if (methods::is(phm, "pheatmap")) {
    grid::grid.draw(methods::slot(phm, "gtable"))
  } else {
    grid::grid.draw(phm)
  }
}


# -----------------------------
# 16. Save outputs
# -----------------------------

write_csv(
  as.data.frame(mat_z) %>%
    rownames_to_column("Protein_ID"),
  "Figure_2A_heatmap_CTRL_EGS_input_matrix.csv"
)

png(
  filename = "Figure_2A_heatmap_CTRL_EGS.png",
  width = 2000,
  height = 2000,
  res = 300,
  type = "cairo-png"
)
grid::grid.newpage()
draw_pheatmap(p_heatmap)
dev.off()

pdf(
  file = "Figure_2A_heatmap_CTRL_EGS.pdf",
  width = 7.5,
  height = 7.5,
  useDingbats = FALSE
)
grid::grid.newpage()
draw_pheatmap(p_heatmap)
dev.off()

svglite::svglite(
  filename = "Figure_2A_heatmap_CTRL_EGS.svg",
  width = 7.5,
  height = 7.5
)
grid::grid.newpage()
draw_pheatmap(p_heatmap)
dev.off()

sink("Figure_2A_heatmap_CTRL_EGS_sessionInfo.txt")
sessionInfo()
sink()

cat("\nSaved:\n")
cat("- Figure_2A_heatmap_CTRL_EGS.png\n")
cat("- Figure_2A_heatmap_CTRL_EGS.pdf\n")
cat("- Figure_2A_heatmap_CTRL_EGS.svg\n")
cat("- Figure_2A_heatmap_CTRL_EGS_input_matrix.csv\n")
cat("- Figure_2A_heatmap_CTRL_EGS_sessionInfo.txt\n")
