# ============================================================
# Figure 2D: PCA of CTRL vs EGS urine proteomes
# ============================================================
#
# Purpose:
# Generate the Figure 2D PCA plot showing sample-level variation
# between CTRL and retained EGS urine samples.
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
#   - Proteins detected in at least 50% of retained samples are kept.
#   - Missing values are imputed using row medians.
#   - Zero-variance proteins are removed.
#   - PCA is performed using centred and scaled data.
#
# Package/method references:
#   - readr::read_tsv:
#     https://readr.tidyverse.org/reference/read_delim.html
#   - dplyr data manipulation:
#     https://dplyr.tidyverse.org/
#   - tibble rowname handling:
#     https://tibble.tidyverse.org/
#   - ggplot2 figure generation:
#     https://ggplot2.tidyverse.org/
#   - ggrepel non-overlapping labels:
#     https://cran.r-project.org/package=ggrepel
#   - svglite SVG export:
#     https://svglite.r-lib.org/
#   - stats::prcomp PCA:
#     https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html
#
# Output files:
#   - Figure_2D_PCA_CTRL_EGS.png
#   - Figure_2D_PCA_CTRL_EGS.pdf
#   - Figure_2D_PCA_CTRL_EGS.svg
#   - Figure_2D_PCA_coordinates.csv
#   - Figure_2D_PCA_sessionInfo.txt
# ============================================================


# ----------------------------
# 0. Load packages
# ----------------------------

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggrepel)
library(svglite)


# ----------------------------
# 1. Read filtered file
# ----------------------------

input_file <- "Filtered_final.tsv"

stopifnot(file.exists(input_file))

Filtered_final <- read_tsv(input_file, show_col_types = FALSE)


# ----------------------------
# 2. Detect PG.Quantity columns
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

cat("\nSample mapping used for Figure 2D:\n")
print(mapping_tbl)

colnames(Filtered_final)[match(sample_cols, colnames(Filtered_final))] <- mapping_tbl$new_name


# ----------------------------
# 4. Define final retained sample order
# ----------------------------

ctrl_keep <- c("CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5")
egs_keep  <- c("EGS1", "EGS2", "EGS3")
keep_samples <- c(ctrl_keep, egs_keep)

stopifnot(all(keep_samples %in% colnames(Filtered_final)))

sample_info <- tibble(
  Sample = keep_samples,
  Group = c(rep("CTRL", length(ctrl_keep)), rep("EGS", length(egs_keep)))
) %>%
  mutate(
    Group = factor(Group, levels = c("CTRL", "EGS"))
  )


# ----------------------------
# 5. Build protein matrix
# ----------------------------

safe_max <- function(x) {
  if (all(is.na(x) | !is.finite(x))) {
    NA_real_
  } else {
    max(x, na.rm = TRUE)
  }
}

mat_tbl <- Filtered_final %>%
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

mat_raw <- mat_tbl %>%
  column_to_rownames("Protein_ID") %>%
  as.matrix()

mode(mat_raw) <- "numeric"


# ----------------------------
# 6. Log2 transform
# ----------------------------

mat_log2 <- log2(mat_raw + 1)


# ----------------------------
# 7. Keep proteins detected in at least 50% of retained samples
# ----------------------------

keep_proteins <- rowSums(!is.na(mat_raw) & is.finite(mat_raw) & mat_raw > 0) >= ceiling(length(keep_samples) / 2)

x_pca <- mat_log2[keep_proteins, keep_samples, drop = FALSE]

cat("\nProteins retained before imputation:", nrow(x_pca), "\n")


# ----------------------------
# 8. Row-median imputation
# ----------------------------

x_imp <- t(apply(x_pca, 1, function(v) {
  v[!is.finite(v) | is.na(v)] <- median(v[is.finite(v)], na.rm = TRUE)
  v
}))

colnames(x_imp) <- colnames(x_pca)
rownames(x_imp) <- rownames(x_pca)


# ----------------------------
# 9. Remove zero-variance proteins
# ----------------------------

x_imp <- x_imp[
  apply(x_imp, 1, sd, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

cat("\nProteins entering PCA:", nrow(x_imp), "\n")


# ----------------------------
# 10. PCA
# ----------------------------

pca <- prcomp(
  t(x_imp),
  center = TRUE,
  scale. = TRUE
)

var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

pc1 <- round(100 * var_expl[1], 1)
pc2 <- round(100 * var_expl[2], 1)

pca_df <- as_tibble(
  pca$x[, 1:2],
  rownames = "Sample"
) %>%
  left_join(sample_info, by = "Sample")

cat("\nPCA variance explained:\n")
cat("PC1:", pc1, "%\n")
cat("PC2:", pc2, "%\n")

print(pca_df)


# ----------------------------
# 11. Colour palette
# ----------------------------

pal_group <- c(
  "CTRL" = "#15a2a2",
  "EGS"  = "#fea802"
)


# ----------------------------
# 12. Plot PCA
# ----------------------------

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, colour = Group)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey70",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey70",
    linewidth = 0.4
  ) +
  geom_point(
    size = 4
  ) +
  geom_text_repel(
    aes(label = Sample),
    size = 3.4,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf,
    show.legend = FALSE,
    seed = 123
  ) +
  scale_colour_manual(
    values = pal_group
  ) +
  labs(
    x = paste0("PC1 (", pc1, "%)"),
    y = paste0("PC2 (", pc2, "%)"),
    colour = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    axis.line = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks = element_line(linewidth = 0.6, colour = "black")
  )

print(p_pca)


# ----------------------------
# 13. Save outputs
# ----------------------------

write_csv(
  pca_df,
  "Figure_2D_PCA_coordinates.csv"
)

ggsave(
  filename = "Figure_2D_PCA_CTRL_EGS.png",
  plot = p_pca,
  width = 7.5,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Figure_2D_PCA_CTRL_EGS.pdf",
  plot = p_pca,
  width = 7.5,
  height = 5.5,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  filename = "Figure_2D_PCA_CTRL_EGS.svg",
  plot = p_pca,
  width = 7.5,
  height = 5.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

sink("Figure_2D_PCA_sessionInfo.txt")
sessionInfo()
sink()


cat("\nSaved:\n")
cat("- Figure_2D_PCA_CTRL_EGS.png\n")
cat("- Figure_2D_PCA_CTRL_EGS.pdf\n")
cat("- Figure_2D_PCA_CTRL_EGS.svg\n")
cat("- Figure_2D_PCA_coordinates.csv\n")
cat("- Figure_2D_PCA_sessionInfo.txt\n")
