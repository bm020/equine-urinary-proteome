# ============================================================
# Figure 2F: Candidate protein boxplots for EGS-increased proteins
# ============================================================
#
# Purpose:
# Generate Figure 2F boxplots showing median-normalised log2
# protein abundance for selected EGS-increased urinary proteins.
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
# Selected proteins:
#   Q95182, ALB, LRG1, VMO1, MSMB, LMAN2
#
# Processing:
#   - Protein quantities are log2-transformed using log2(x + 1).
#   - Sample-wise median normalisation is applied.
#   - Missing PG.Genes values are replaced with Protein_ID.
#   - Duplicate gene entries are collapsed per sample using maximum abundance.
#   - Boxplots are generated without individual data points.
#
# Package/method references:
#   - readr::read_tsv:
#     https://readr.tidyverse.org/reference/read_delim.html
#   - dplyr data manipulation:
#     https://dplyr.tidyverse.org/
#   - tidyr data reshaping:
#     https://tidyr.tidyverse.org/
#   - tibble table handling:
#     https://tibble.tidyverse.org/
#   - stringr string handling:
#     https://stringr.tidyverse.org/
#   - Base R boxplot:
#     https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/boxplot.html
#   - grDevices SVG export:
#     https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/svg.html
#
# Output files:
#   - Figure_2F_candidate_boxplots.png
#   - Figure_2F_candidate_boxplots.pdf
#   - Figure_2F_candidate_boxplots.svg
#   - Figure_2F_candidate_boxplots_selected_proteins.csv
#   - Figure_2F_candidate_boxplots_sessionInfo.txt
# ============================================================


# ----------------------------
# 0. Load packages
# ----------------------------

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)


# ----------------------------
# 1. Read filtered file
# ----------------------------

input_file <- "Filtered_final.tsv"

stopifnot(file.exists(input_file))

Filtered_final <- read_tsv(input_file, show_col_types = FALSE)


# ----------------------------
# 2. Detect PG.Quantity columns
# ----------------------------

quantity_cols <- grep(
  "Horse_urine_U.*PG\\.Quantity$",
  colnames(Filtered_final),
  value = TRUE
)

stopifnot(length(quantity_cols) == 10)

cat("\nDetected PG.Quantity columns:\n")
print(quantity_cols)


# ----------------------------
# 3. Build working table and rename samples
# ----------------------------

data_wide <- Filtered_final %>%
  select(
    Protein_ID = PG.ProteinGroups,
    Gene = PG.Genes,
    all_of(quantity_cols)
  )

colnames(data_wide) <- case_when(
  colnames(data_wide) == "Protein_ID" ~ "Protein_ID",
  colnames(data_wide) == "Gene" ~ "Gene",

  str_detect(colnames(data_wide), "Horse_urine_U10_") ~ "EGS5_excluded",
  str_detect(colnames(data_wide), "Horse_urine_U9_")  ~ "EGS4_excluded",
  str_detect(colnames(data_wide), "Horse_urine_U8_")  ~ "EGS3",
  str_detect(colnames(data_wide), "Horse_urine_U7_")  ~ "EGS2",
  str_detect(colnames(data_wide), "Horse_urine_U6_")  ~ "EGS1",

  str_detect(colnames(data_wide), "Horse_urine_U5_")  ~ "CTRL1",
  str_detect(colnames(data_wide), "Horse_urine_U4_")  ~ "CTRL2",
  str_detect(colnames(data_wide), "Horse_urine_U3_")  ~ "CTRL3",
  str_detect(colnames(data_wide), "Horse_urine_U2_")  ~ "CTRL4",
  str_detect(colnames(data_wide), "Horse_urine_U1_")  ~ "CTRL5",

  TRUE ~ colnames(data_wide)
)


# ----------------------------
# 4. Define retained sample order
# ----------------------------

ctrl_keep <- c("CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5")
egs_keep  <- c("EGS1", "EGS2", "EGS3")
samples_keep <- c(ctrl_keep, egs_keep)

stopifnot(all(samples_keep %in% colnames(data_wide)))

cat("\nSamples retained for Figure 2F boxplots:\n")
print(samples_keep)

cat("\nDetected protein groups per retained sample in Filtered_final.tsv:\n")
print(
  colSums(
    !is.na(data_wide[, samples_keep]) &
      is.finite(as.matrix(data_wide[, samples_keep])) &
      as.matrix(data_wide[, samples_keep]) > 0
  )
)

data_wide <- data_wide %>%
  select(
    Protein_ID,
    Gene,
    all_of(samples_keep)
  )


# ----------------------------
# 5. Clean gene symbols
# ----------------------------
# If PG.Genes is missing, fall back to Protein_ID.
# This preserves accession-style entries such as Q95182.
# ----------------------------

data_wide <- data_wide %>%
  mutate(
    Protein_ID = as.character(Protein_ID),
    Gene = as.character(Gene),
    Gene = if_else(is.na(Gene) | Gene == "", Protein_ID, Gene),
    Gene = str_replace_all(Gene, ";", ",")
  ) %>%
  separate_rows(Gene, sep = ",") %>%
  mutate(
    Gene = str_trim(Gene)
  ) %>%
  filter(
    !is.na(Gene),
    Gene != ""
  )


# ----------------------------
# 6. Convert to long format and log2-transform
# ----------------------------

data_long <- data_wide %>%
  pivot_longer(
    cols = all_of(samples_keep),
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
# 7. Median-normalise within each sample
# ----------------------------

data_long <- data_long %>%
  group_by(Sample) %>%
  mutate(
    median_normalised_log2 =
      log2_abundance - median(log2_abundance, na.rm = TRUE)
  ) %>%
  ungroup()


# ----------------------------
# 8. Collapse duplicate gene entries per sample
# ----------------------------

safe_max <- function(x) {
  if (all(is.na(x) | !is.finite(x))) {
    NA_real_
  } else {
    max(x, na.rm = TRUE)
  }
}

data_gene <- data_long %>%
  group_by(Gene, Sample, Group) %>%
  summarise(
    abundance = safe_max(median_normalised_log2),
    .groups = "drop"
  ) %>%
  filter(is.finite(abundance))


# ----------------------------
# 9. Define selected candidate proteins
# ----------------------------

selected_genes <- c(
  "Q95182",
  "ALB",
  "LRG1",
  "VMO1",
  "MSMB",
  "LMAN2"
)

display_label_map <- c(
  "Q95182" = "Q95182",
  "ALB"   = "ALB",
  "LRG1"  = "LRG1",
  "VMO1"  = "VMO1",
  "MSMB"  = "MSMB",
  "LMAN2" = "LMAN2"
)

missing_genes <- setdiff(selected_genes, unique(data_gene$Gene))

if (length(missing_genes) > 0) {
  warning(
    "The following selected proteins were not found in data_gene: ",
    paste(missing_genes, collapse = ", ")
  )
}

selected_genes_found <- selected_genes[selected_genes %in% unique(data_gene$Gene)]

if (length(selected_genes_found) == 0) {
  stop("None of the selected proteins were found in the data.")
}


# ----------------------------
# 10. Calculate summary fold-change table
# ----------------------------

fc_table <- data_gene %>%
  filter(Gene %in% selected_genes_found) %>%
  group_by(Gene, Group) %>%
  summarise(
    median_abundance = median(abundance, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(median_abundance, n)
  ) %>%
  mutate(
    log2FC_EGS_vs_CTRL = median_abundance_EGS - median_abundance_CTRL,
    FC_EGS_vs_CTRL = 2^log2FC_EGS_vs_CTRL,
    Display_label = unname(display_label_map[Gene])
  ) %>%
  arrange(match(Gene, selected_genes_found))

cat("\nSelected proteins for Figure 2F boxplots:\n")
print(fc_table)

write_csv(
  fc_table,
  "Figure_2F_candidate_boxplots_selected_proteins.csv"
)


# ----------------------------
# 11. Prepare plotting data
# ----------------------------

plot_data <- data_gene %>%
  filter(Gene %in% selected_genes_found) %>%
  mutate(
    Gene = factor(Gene, levels = selected_genes_found),
    Display_label = unname(display_label_map[as.character(Gene)]),
    Display_label = factor(
      Display_label,
      levels = unname(display_label_map[selected_genes_found])
    ),
    Group = factor(Group, levels = c("CTRL", "EGS"))
  )


# ----------------------------
# 12. Define colours
# ----------------------------

group_cols <- c(
  "CTRL" = "#15a2a2",
  "EGS"  = "#fea802"
)


# ----------------------------
# 13. Plotting function
# ----------------------------

plot_boxplot <- function() {

  # Two panels:
  # left = boxplots
  # right = legend outside the plot area

  layout(
    matrix(c(1, 2), nrow = 1),
    widths = c(5.0, 1.2)
  )

  # Main boxplot panel
  par(
    mar = c(6.2, 5.5, 1.2, 0.5),
    family = "sans"
  )

  bp <- boxplot(
    abundance ~ Group * Display_label,
    data = plot_data,
    plot = FALSE,
    outline = FALSE
  )

  y_range <- range(plot_data$abundance, na.rm = TRUE)
  y_padding <- 0.10 * diff(y_range)

  if (y_padding == 0 || is.na(y_padding)) {
    y_padding <- 1
  }

  y_lim <- c(
    y_range[1] - y_padding,
    y_range[2] + y_padding
  )

  plot(
    NA,
    xlim = c(0.5, length(bp$names) + 0.5),
    ylim = y_lim,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "o"
  )

  # Grey separator lines between protein groups
  for (i in seq(0.5, length(bp$names) + 0.5, by = 2)) {
    abline(
      v = i,
      lty = 1,
      col = "grey85",
      lwd = 1
    )
  }

  boxplot(
    abundance ~ Group * Display_label,
    data = plot_data,
    boxwex = 0.45,
    col = rep(group_cols, times = length(selected_genes_found)),
    xaxt = "n",
    yaxt = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    outline = FALSE,
    border = "black",
    lwd = 1.4,
    add = TRUE
  )

  axis(
    side = 2,
    las = 1,
    cex.axis = 1.1,
    lwd = 1.2,
    lwd.ticks = 1.2,
    col.axis = "black"
  )

  axis(
    side = 1,
    at = seq(1.5, length(bp$names), by = 2),
    labels = unname(display_label_map[selected_genes_found]),
    tick = FALSE,
    cex.axis = 1.05,
    las = 1,
    col.axis = "black"
  )

  box(
    lwd = 1.2,
    col = "black"
  )

  title(
    ylab = expression("Median-normalised log"[2] * " protein abundance"),
    cex.lab = 1.25,
    line = 3.6
  )

  # Legend panel
  par(
    mar = c(6.2, 0.5, 1.2, 1.5)
  )

  plot.new()

  legend(
    "center",
    legend = c("CTRL", "EGS"),
    fill = group_cols,
    border = "black",
    bty = "n",
    cex = 1.2
  )
}


# ----------------------------
# 14. Export figures
# ----------------------------

png(
  filename = "Figure_2F_candidate_boxplots.png",
  width = 9.5,
  height = 6,
  units = "in",
  res = 600,
  type = "cairo-png"
)
plot_boxplot()
dev.off()

pdf(
  file = "Figure_2F_candidate_boxplots.pdf",
  width = 9.5,
  height = 6,
  useDingbats = FALSE
)
plot_boxplot()
dev.off()

svg(
  filename = "Figure_2F_candidate_boxplots.svg",
  width = 9.5,
  height = 6
)
plot_boxplot()
dev.off()


# ----------------------------
# 15. Save session information
# ----------------------------

sink("Figure_2F_candidate_boxplots_sessionInfo.txt")
sessionInfo()
sink()


cat("\nSaved:\n")
cat("- Figure_2F_candidate_boxplots.png\n")
cat("- Figure_2F_candidate_boxplots.pdf\n")
cat("- Figure_2F_candidate_boxplots.svg\n")
cat("- Figure_2F_candidate_boxplots_selected_proteins.csv\n")
cat("- Figure_2F_candidate_boxplots_sessionInfo.txt\n")
