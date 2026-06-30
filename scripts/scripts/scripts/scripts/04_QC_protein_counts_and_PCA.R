# ============================================================
# Supplementary Figure S1: QC protein counts and PCA
# ============================================================
#
# Manuscript:
# Defining the Equine Urinary Proteome: A Reference Baseline
# for Biomarker Discovery
#
# Purpose:
# This script generates quality-control figures summarising:
#   1. Protein identifications per sample.
#   2. Protein identifications by sample group.
#   3. PCA of urine proteomics profiles across fresh, CTRL and EGS samples.
#
# Input files:
#   - Horse_urine_U1&U3_horse_fasta_020625_report.pg_matrix.tsv
#   - Horse_urine_U2_donkey_fasta_160525_report.pg_matrix.tsv
#   - Filtered_final.tsv
#
# Notes:
#   - Fresh urine samples: U1, U2 and U3.
#   - Archived control samples: CTRL1-CTRL5.
#   - Archived EGS samples: EGS1-EGS5.
#   - EGS4 and EGS5 are included here only for QC visualisation.
#     They were excluded from the final downstream comparative analysis
#     due to low proteome coverage.
#
# Package/method references:
#   - readr::read_tsv for importing TSV files:
#     https://readr.tidyverse.org/reference/read_delim.html
#   - dplyr::group_by and dplyr::summarise:
#     https://dplyr.tidyverse.org/reference/group_by.html
#     https://dplyr.tidyverse.org/reference/summarise.html
#   - ggplot2 for figure generation:
#     https://ggplot2.tidyverse.org/
#   - ggsave for exporting figures:
#     https://ggplot2.tidyverse.org/reference/ggsave.html
#   - stats::prcomp for PCA:
#     https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html
# ============================================================


# ----------------------------
# 0. Load packages
# ----------------------------

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(grid)


# ============================================================
# Part A: Protein identifications per sample
# ============================================================

# ----------------------------
# 1. QC protein counts
# ----------------------------

qc_counts <- tibble::tribble(
  ~Sample,  ~Group,   ~Protein_IDs,
  "U1",     "Fresh",  1094,
  "U2",     "Fresh",  1006,
  "U3",     "Fresh",  1136,
  "CTRL1",  "CTRL",    556,
  "CTRL2",  "CTRL",    538,
  "CTRL3",  "CTRL",    631,
  "CTRL4",  "CTRL",    520,
  "CTRL5",  "CTRL",    661,
  "EGS1",   "EGS",     678,
  "EGS2",   "EGS",     732,
  "EGS3",   "EGS",     613,
  "EGS4",   "EGS",      81,
  "EGS5",   "EGS",     255
)


# ----------------------------
# 2. Order samples and add group spacing
# ----------------------------

qc_sample_plot_df <- qc_counts %>%
  mutate(
    Group = factor(Group, levels = c("Fresh", "CTRL", "EGS")),
    Sample = factor(
      Sample,
      levels = c(
        "U1", "U2", "U3",
        "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5",
        "EGS1", "EGS2", "EGS3", "EGS4", "EGS5"
      )
    ),
    xpos = case_when(
      Sample == "U1"    ~ 1,
      Sample == "U2"    ~ 2,
      Sample == "U3"    ~ 3,
      Sample == "CTRL1" ~ 5,
      Sample == "CTRL2" ~ 6,
      Sample == "CTRL3" ~ 7,
      Sample == "CTRL4" ~ 8,
      Sample == "CTRL5" ~ 9,
      Sample == "EGS1"  ~ 11,
      Sample == "EGS2"  ~ 12,
      Sample == "EGS3"  ~ 13,
      Sample == "EGS4"  ~ 14,
      Sample == "EGS5"  ~ 15
    )
  )


# ----------------------------
# 3. Colours
# ----------------------------

group_cols <- c(
  "Fresh" = "#0450B4",
  "CTRL"  = "#15A2A2",
  "EGS"   = "#FEA802"
)


# ----------------------------
# 4. Plot protein identifications per sample
# ----------------------------

p_counts_per_sample <- ggplot(
  qc_sample_plot_df,
  aes(x = xpos, y = Protein_IDs, fill = Group)
) +
  geom_col(
    width = 0.72,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = Protein_IDs),
    vjust = -0.35,
    size = 3.8,
    colour = "black"
  ) +
  scale_fill_manual(values = group_cols, breaks = c("Fresh", "CTRL", "EGS")) +
  scale_x_continuous(
    breaks = qc_sample_plot_df$xpos,
    labels = qc_sample_plot_df$Sample,
    expand = expansion(mult = c(0.015, 0.08))
  ) +
  scale_y_continuous(
    limits = c(0, 1300),
    breaks = c(0, 250, 500, 750, 1000, 1250),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Proteins detected (n)",
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      size = 10.5,
      colour = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.title.y = element_text(
      size = 13,
      colour = "black",
      margin = margin(r = 9)
    ),
    legend.position = "right",
    legend.text = element_text(size = 11, colour = "black"),
    legend.key.size = unit(0.7, "cm"),
    axis.line = element_line(linewidth = 0.55, colour = "black"),
    axis.ticks = element_line(linewidth = 0.55, colour = "black"),
    plot.margin = margin(10, 20, 15, 10)
  )

print(p_counts_per_sample)


# ----------------------------
# 5. Save protein identifications per sample
# ----------------------------

ggsave(
  "Supplementary_Figure_S1A_protein_detection_per_sample.png",
  plot = p_counts_per_sample,
  width = 10.5,
  height = 5.6,
  units = "in",
  dpi = 600
)

ggsave(
  "Supplementary_Figure_S1A_protein_detection_per_sample.svg",
  plot = p_counts_per_sample,
  width = 10.5,
  height = 5.6,
  units = "in"
)

ggsave(
  "Supplementary_Figure_S1A_protein_detection_per_sample.pdf",
  plot = p_counts_per_sample,
  width = 10.5,
  height = 5.6,
  units = "in",
  useDingbats = FALSE
)


# ============================================================
# Part B: Protein identifications by group
# ============================================================

# ----------------------------
# 6. Prepare data for group-level labels
# ----------------------------

qc_group_plot_df <- qc_counts %>%
  mutate(
    Group = factor(Group, levels = c("Fresh", "CTRL", "EGS")),
    group_num = as.numeric(Group),
    label_x = case_when(
      Group == "Fresh" ~ group_num + 0.32,
      Group == "CTRL"  ~ group_num - 0.45,
      Group == "EGS"   ~ group_num + 0.32
    ),
    label_y = case_when(
      Sample == "CTRL5" ~ Protein_IDs + 35,
      Sample == "CTRL3" ~ Protein_IDs + 18,
      Sample == "CTRL1" ~ Protein_IDs + 5,
      Sample == "CTRL2" ~ Protein_IDs - 12,
      Sample == "CTRL4" ~ Protein_IDs - 30,
      TRUE ~ Protein_IDs
    ),
    label_hjust = case_when(
      Group == "CTRL" ~ 1,
      TRUE ~ 0
    )
  )


# ----------------------------
# 7. Plot protein identifications by group
# ----------------------------

p_counts_by_group <- ggplot(
  qc_group_plot_df,
  aes(x = Group, y = Protein_IDs, fill = Group)
) +
  geom_boxplot(
    width = 0.28,
    alpha = 0.25,
    colour = "black",
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  geom_point(
    aes(colour = Group),
    size = 3,
    stroke = 0.7
  ) +
  geom_segment(
    aes(
      x = group_num,
      xend = label_x,
      y = Protein_IDs,
      yend = label_y,
      colour = Group
    ),
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_text(
    aes(
      x = label_x,
      y = label_y,
      label = Sample,
      colour = Group,
      hjust = label_hjust
    ),
    size = 3.4,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_cols) +
  scale_colour_manual(values = group_cols) +
  scale_y_continuous(
    limits = c(0, 1250),
    breaks = seq(0, 1250, 250),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Protein identifications (n)"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    plot.margin = margin(10, 90, 10, 90)
  )

print(p_counts_by_group)


# ----------------------------
# 8. Save protein identifications by group
# ----------------------------

ggsave(
  "Supplementary_Figure_S1B_protein_identifications_by_group.svg",
  plot = p_counts_by_group,
  width = 7.5,
  height = 4.2,
  units = "in"
)

ggsave(
  "Supplementary_Figure_S1B_protein_identifications_by_group.pdf",
  plot = p_counts_by_group,
  width = 7.5,
  height = 4.2,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  "Supplementary_Figure_S1B_protein_identifications_by_group.png",
  plot = p_counts_by_group,
  width = 7.5,
  height = 4.2,
  units = "in",
  dpi = 600
)


# ============================================================
# Part C: QC PCA
# ============================================================

# ----------------------------
# 9. File paths
# ----------------------------

u1u3_file  <- "Horse_urine_U1&U3_horse_fasta_020625_report.pg_matrix.tsv"
u2_file    <- "Horse_urine_U2_donkey_fasta_160525_report.pg_matrix.tsv"
final_file <- "Filtered_final.tsv"

stopifnot(file.exists(u1u3_file))
stopifnot(file.exists(u2_file))
stopifnot(file.exists(final_file))


# ----------------------------
# 10. Read files
# ----------------------------

U1U3  <- read_tsv(u1u3_file, show_col_types = FALSE)
U2    <- read_tsv(u2_file, show_col_types = FALSE)
FINAL <- read_tsv(final_file, show_col_types = FALSE)


# ----------------------------
# 11. Sample metadata
# ----------------------------

sample_info <- tibble::tribble(
  ~Sample,  ~Group,
  "U1",     "Fresh",
  "U2",     "Fresh",
  "U3",     "Fresh",
  "CTRL1",  "CTRL",
  "CTRL2",  "CTRL",
  "CTRL3",  "CTRL",
  "CTRL4",  "CTRL",
  "CTRL5",  "CTRL",
  "EGS1",   "EGS",
  "EGS2",   "EGS",
  "EGS3",   "EGS",
  "EGS4",   "EGS",
  "EGS5",   "EGS"
) %>%
  mutate(
    Group = factor(Group, levels = c("Fresh", "CTRL", "EGS")),
    Sample = factor(
      Sample,
      levels = c(
        "U1", "U2", "U3",
        "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5",
        "EGS1", "EGS2", "EGS3", "EGS4", "EGS5"
      )
    )
  )


# ----------------------------
# 12. Extract U1 and U3 data
# ----------------------------

u1u3_intensity_cols <- grep("\\.wiff$", colnames(U1U3), value = TRUE)

u1u3_intensity_cols <- u1u3_intensity_cols[
  !str_detect(u1u3_intensity_cols, regex("pooled", ignore_case = TRUE))
]

U1U3_mat <- U1U3 %>%
  select(
    Protein_ID = Protein.Group,
    all_of(u1u3_intensity_cols)
  )

colnames(U1U3_mat) <- case_when(
  colnames(U1U3_mat) == "Protein_ID" ~ "Protein_ID",
  str_detect(colnames(U1U3_mat), "Horse_urine_U1_") ~ "U1",
  str_detect(colnames(U1U3_mat), "Horse_urine_U3_") ~ "U3",
  TRUE ~ colnames(U1U3_mat)
)


# ----------------------------
# 13. Extract U2 data
# ----------------------------

u2_intensity_cols <- grep("\\.wiff$", colnames(U2), value = TRUE)

U2_mat <- U2 %>%
  select(
    Protein_ID = Protein.Group,
    all_of(u2_intensity_cols)
  )

colnames(U2_mat) <- case_when(
  colnames(U2_mat) == "Protein_ID" ~ "Protein_ID",
  str_detect(colnames(U2_mat), "Horse_urine_U2_") ~ "U2",
  TRUE ~ colnames(U2_mat)
)


# ----------------------------
# 14. Extract CTRL and EGS data
# ----------------------------

final_quantity_cols <- grep("PG\\.Quantity$", colnames(FINAL), value = TRUE)

FINAL_mat <- FINAL %>%
  select(
    Protein_ID = PG.ProteinGroups,
    all_of(final_quantity_cols)
  )

# Mapping:
# U10 = EGS1
# U9  = EGS2
# U8  = EGS3
# U7  = EGS4
# U6  = EGS5
# U5  = CTRL1
# U4  = CTRL2
# U3  = CTRL3
# U2  = CTRL4
# U1  = CTRL5

colnames(FINAL_mat) <- case_when(
  colnames(FINAL_mat) == "Protein_ID" ~ "Protein_ID",
  str_detect(colnames(FINAL_mat), "Horse_urine_U10_") ~ "EGS1",
  str_detect(colnames(FINAL_mat), "Horse_urine_U9_")  ~ "EGS2",
  str_detect(colnames(FINAL_mat), "Horse_urine_U8_")  ~ "EGS3",
  str_detect(colnames(FINAL_mat), "Horse_urine_U7_")  ~ "EGS4",
  str_detect(colnames(FINAL_mat), "Horse_urine_U6_")  ~ "EGS5",
  str_detect(colnames(FINAL_mat), "Horse_urine_U5_")  ~ "CTRL1",
  str_detect(colnames(FINAL_mat), "Horse_urine_U4_")  ~ "CTRL2",
  str_detect(colnames(FINAL_mat), "Horse_urine_U3_")  ~ "CTRL3",
  str_detect(colnames(FINAL_mat), "Horse_urine_U2_")  ~ "CTRL4",
  str_detect(colnames(FINAL_mat), "Horse_urine_U1_")  ~ "CTRL5",
  TRUE ~ colnames(FINAL_mat)
)


# ----------------------------
# 15. Check extracted columns
# ----------------------------

cat("\nU1/U3 columns:\n")
print(colnames(U1U3_mat))

cat("\nU2 columns:\n")
print(colnames(U2_mat))

cat("\nFiltered_final columns:\n")
print(colnames(FINAL_mat))


# ----------------------------
# 16. Merge by protein ID
# ----------------------------

combined_df <- U1U3_mat %>%
  full_join(U2_mat, by = "Protein_ID") %>%
  full_join(FINAL_mat, by = "Protein_ID")


# Collapse duplicated protein IDs if present
combined_df <- combined_df %>%
  group_by(Protein_ID) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    across(where(is.numeric), ~ ifelse(is.nan(.), NA, .))
  )


# ----------------------------
# 17. Keep and order samples
# ----------------------------

expected_samples <- c(
  "U1", "U2", "U3",
  "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5",
  "EGS1", "EGS2", "EGS3", "EGS4", "EGS5"
)

combined_df <- combined_df %>%
  select(Protein_ID, any_of(expected_samples))

cat("\nSamples found for PCA:\n")
print(colnames(combined_df))


# ----------------------------
# 18. Prepare intensity matrix
# ----------------------------

combined_matrix <- combined_df %>%
  select(-Protein_ID) %>%
  mutate(across(everything(), ~ as.numeric(.)))


# ----------------------------
# 19. Log2 transform if needed
# ----------------------------

max_value <- max(as.matrix(combined_matrix), na.rm = TRUE)

if (max_value > 100) {
  combined_matrix <- log2(combined_matrix + 1)
  cat("\nData were log2-transformed using log2(x + 1).\n")
} else {
  cat("\nData already appear to be log-scaled; no log2 transformation applied.\n")
}


# ----------------------------
# 20. Filter proteins by completeness
# ----------------------------

keep_proteins <- rowMeans(!is.na(combined_matrix)) >= 0.50
pca_matrix <- combined_matrix[keep_proteins, ]

cat("\nProteins retained for PCA:", nrow(pca_matrix), "\n")


# ----------------------------
# 21. Impute missing values by row median
# ----------------------------

pca_matrix_imputed <- pca_matrix %>%
  as.data.frame() %>%
  mutate(row_median = apply(., 1, function(x) median(x, na.rm = TRUE))) %>%
  mutate(
    across(
      -row_median,
      ~ ifelse(is.na(.), row_median, .)
    )
  ) %>%
  select(-row_median)


# ----------------------------
# 22. Median-normalise per sample
# ----------------------------

pca_matrix_norm <- sweep(
  pca_matrix_imputed,
  2,
  apply(pca_matrix_imputed, 2, median, na.rm = TRUE),
  FUN = "-"
)


# ----------------------------
# 23. Run PCA
# ----------------------------

pca_input <- t(pca_matrix_norm)

cat("\nPCA input dimensions:\n")
cat("Samples:", nrow(pca_input), "\n")
cat("Proteins:", ncol(pca_input), "\n")

pca_res <- prcomp(
  pca_input,
  center = TRUE,
  scale. = TRUE
)

pca_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("Sample") %>%
  left_join(sample_info, by = "Sample") %>%
  filter(!is.na(Group))

print(pca_df)


# ----------------------------
# 24. Manual label positions
# ----------------------------

label_df <- pca_df %>%
  mutate(
    dx = case_when(
      Sample == "U1"    ~  2.4,
      Sample == "U3"    ~ -2.0,
      Sample == "U2"    ~  1.3,
      Sample == "CTRL1" ~ -1.2,
      Sample == "CTRL2" ~  1.3,
      Sample == "CTRL3" ~ -1.8,
      Sample == "CTRL4" ~  1.0,
      Sample == "CTRL5" ~ -1.0,
      Sample == "EGS1"  ~  1.5,
      Sample == "EGS2"  ~ -1.3,
      Sample == "EGS3"  ~ -0.6,
      Sample == "EGS4"  ~ -0.3,
      Sample == "EGS5"  ~ -0.5,
      TRUE ~ 0
    ),
    dy = case_when(
      Sample == "U1"    ~ -0.3,
      Sample == "U3"    ~  1.2,
      Sample == "U2"    ~ -0.4,
      Sample == "CTRL1" ~  1.3,
      Sample == "CTRL2" ~ -2.2,
      Sample == "CTRL3" ~ -1.5,
      Sample == "CTRL4" ~ -0.7,
      Sample == "CTRL5" ~  1.0,
      Sample == "EGS1"  ~  1.4,
      Sample == "EGS2"  ~ -1.9,
      Sample == "EGS3"  ~  6.4,
      Sample == "EGS4"  ~ -2.7,
      Sample == "EGS5"  ~ -1.7,
      TRUE ~ 0
    ),
    label_x = PC1 + dx,
    label_y = PC2 + dy,
    hjust = case_when(
      dx > 0 ~ 0,
      dx < 0 ~ 1,
      TRUE ~ 0.5
    )
  )


# ----------------------------
# 25. Plot PCA
# ----------------------------

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    colour = "grey75",
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3,
    colour = "grey75",
    linetype = "dashed"
  ) +
  geom_point(
    aes(fill = Group),
    shape = 21,
    size = 4.2,
    stroke = 0.45,
    colour = "black",
    alpha = 0.95
  ) +
  geom_segment(
    data = label_df,
    aes(
      x = PC1,
      y = PC2,
      xend = label_x,
      yend = label_y,
      colour = Group
    ),
    linewidth = 0.35,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_text(
    data = label_df,
    aes(
      x = label_x,
      y = label_y,
      label = Sample,
      colour = Group,
      hjust = hjust
    ),
    size = 4.4,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = group_cols,
    breaks = c("Fresh", "CTRL", "EGS")
  ) +
  scale_colour_manual(
    values = group_cols,
    breaks = c("Fresh", "CTRL", "EGS")
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.08, 0.16))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.10, 0.12))
  ) +
  labs(
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 16, colour = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 13, colour = "black"),
    legend.key.size = unit(0.75, "cm"),
    axis.line = element_line(linewidth = 0.65, colour = "black"),
    axis.ticks = element_line(linewidth = 0.65, colour = "black"),
    plot.margin = margin(10, 28, 12, 10)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, size = 4.5, colour = "black")
    )
  )

print(p_pca)


# ----------------------------
# 26. Save PCA
# ----------------------------

ggsave(
  "Supplementary_Figure_S1C_QC_PCA_urine_protein_intensity.png",
  plot = p_pca,
  width = 7.0,
  height = 5.2,
  units = "in",
  dpi = 600
)

ggsave(
  "Supplementary_Figure_S1C_QC_PCA_urine_protein_intensity.svg",
  plot = p_pca,
  width = 7.0,
  height = 5.2,
  units = "in"
)

ggsave(
  "Supplementary_Figure_S1C_QC_PCA_urine_protein_intensity.pdf",
  plot = p_pca,
  width = 7.0,
  height = 5.2,
  units = "in",
  useDingbats = FALSE
)


# ----------------------------
# 27. Final message
# ----------------------------

cat("\nSaved:\n")
cat("- Supplementary_Figure_S1A_protein_detection_per_sample.png/.svg/.pdf\n")
cat("- Supplementary_Figure_S1B_protein_identifications_by_group.png/.svg/.pdf\n")
cat("- Supplementary_Figure_S1C_QC_PCA_urine_protein_intensity.png/.svg/.pdf\n")
