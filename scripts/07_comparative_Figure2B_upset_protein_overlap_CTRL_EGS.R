# ============================================================
# Figure 2B: UpSet plot of protein detection overlap in CTRL vs EGS urine
# ============================================================
#
# Purpose:
# Generate the Figure 2B UpSet plot showing overlap of identified
# protein groups across individual CTRL and retained EGS urine samples.
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
#   - UpSetR set intersection visualisation:
#     https://cran.r-project.org/package=UpSetR
#
# Output files:
#   - Figure_2B_upset_protein_overlap_CTRL_EGS.png
#   - Figure_2B_upset_protein_overlap_CTRL_EGS.pdf
#   - Figure_2B_upset_protein_overlap_CTRL_EGS.svg
#   - Figure_2B_upset_protein_counts.csv
#   - Figure_2B_upset_sessionInfo.txt
# ============================================================


# ----------------------------
# 0. Load packages
# ----------------------------

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(UpSetR)
library(svglite)


# ----------------------------
# 1. Read filtered file
# ----------------------------

input_file <- "Filtered_final.tsv"

stopifnot(file.exists(input_file))

Filtered_final <- read_tsv(input_file, show_col_types = FALSE)


# ----------------------------
# 2. Detect PG.Quantity sample columns
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
# 3. Rename columns using correct biological mapping
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

cat("\nSample mapping used for Figure 2B:\n")
print(mapping_tbl)

colnames(Filtered_final)[match(sample_cols, colnames(Filtered_final))] <- mapping_tbl$new_name


# ----------------------------
# 4. Define final retained sample order
# ----------------------------

sample_order <- c(
  "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5",
  "EGS1", "EGS2", "EGS3"
)

stopifnot(all(sample_order %in% colnames(Filtered_final)))


# ----------------------------
# 5. Define colours
# ----------------------------

pal_group <- c(
  "CTRL" = "#15a2a2",
  "EGS"  = "#fea802"
)

sets_bar_cols <- ifelse(
  grepl("^CTRL", sample_order),
  pal_group["CTRL"],
  pal_group["EGS"]
)


# ----------------------------
# 6. Build protein detection lists
# ----------------------------

protein_ids <- Filtered_final$PG.ProteinGroups

if (is.null(protein_ids)) {
  stop("Column PG.ProteinGroups was not found in Filtered_final.tsv.")
}

x <- Filtered_final[, sample_order, drop = FALSE]

pres <- as.data.frame(
  lapply(x, function(v) {
    v <- as.numeric(v)
    !is.na(v) & is.finite(v) & v > 0
  })
)

prot_list <- lapply(sample_order, function(s) {
  unique(protein_ids[pres[[s]]])
})

names(prot_list) <- sample_order

upset_input <- UpSetR::fromList(prot_list)


# ----------------------------
# 7. Save protein count table
# ----------------------------

protein_counts <- tibble(
  Sample = sample_order,
  Group = if_else(str_detect(sample_order, "^CTRL"), "CTRL", "EGS"),
  Detected_protein_groups = sapply(prot_list, length)
)

cat("\nDetected protein groups per retained sample:\n")
print(protein_counts)

write_csv(
  protein_counts,
  "Figure_2B_upset_protein_counts.csv"
)


# ----------------------------
# 8. Plot function
# ----------------------------

plot_upset <- function() {
  UpSetR::upset(
    upset_input,
    sets = sample_order,
    nsets = length(sample_order),
    keep.order = TRUE,
    order.by = "freq",
    decreasing = TRUE,
    nintersects = 30,
    mb.ratio = c(0.65, 0.35),
    number.angles = 0,
    text.scale = c(1.7, 1.4, 1.2, 1.2, 1.2, 1.2),
    point.size = 3.1,
    line.size = 1.1,
    main.bar.color = "grey30",
    sets.bar.color = sets_bar_cols,
    matrix.color = "grey30",
    mainbar.y.label = "Intersection size",
    sets.x.label = "Detected proteins"
  )
}


# ----------------------------
# 9. Save figure outputs
# ----------------------------

png(
  filename = "Figure_2B_upset_protein_overlap_CTRL_EGS.png",
  width = 3200,
  height = 1700,
  res = 300,
  type = "cairo-png"
)
plot_upset()
dev.off()

pdf(
  file = "Figure_2B_upset_protein_overlap_CTRL_EGS.pdf",
  width = 11,
  height = 6.5,
  useDingbats = FALSE
)
plot_upset()
dev.off()

svglite::svglite(
  filename = "Figure_2B_upset_protein_overlap_CTRL_EGS.svg",
  width = 11,
  height = 6.5
)
plot_upset()
dev.off()


# ----------------------------
# 10. Save session information
# ----------------------------

sink("Figure_2B_upset_sessionInfo.txt")
sessionInfo()
sink()


cat("\nSaved:\n")
cat("- Figure_2B_upset_protein_overlap_CTRL_EGS.png\n")
cat("- Figure_2B_upset_protein_overlap_CTRL_EGS.pdf\n")
cat("- Figure_2B_upset_protein_overlap_CTRL_EGS.svg\n")
cat("- Figure_2B_upset_protein_counts.csv\n")
cat("- Figure_2B_upset_sessionInfo.txt\n")
