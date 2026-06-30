# ============================================================
# Figure 2E: Volcano plot of CTRL vs EGS urine proteomes
# ============================================================
#
# Purpose:
# Generate the Figure 2E volcano plot showing proteins with
# higher abundance in CTRL or EGS urine samples.
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
#   - Proteins require at least 2 valid values in each group.
#   - Zero-variance proteins are removed.
#   - Differential abundance is assessed using limma.
#   - Thresholds used for plotting:
#       log2FC >= 0.8 and p-value < 0.1: higher in EGS
#       log2FC <= -0.8 and p-value < 0.1: higher in CTRL
#
# Package/method references:
#   - readr::read_tsv:
#     https://readr.tidyverse.org/reference/read_delim.html
#   - dplyr data manipulation:
#     https://dplyr.tidyverse.org/
#   - tibble rowname handling:
#     https://tibble.tidyverse.org/
#   - stringr string handling:
#     https://stringr.tidyverse.org/
#   - limma empirical Bayes linear modelling:
#     https://bioconductor.org/packages/release/bioc/html/limma.html
#   - ggplot2 figure generation:
#     https://ggplot2.tidyverse.org/
#   - ggrepel non-overlapping labels:
#     https://cran.r-project.org/package=ggrepel
#   - svglite SVG export:
#     https://svglite.r-lib.org/
#
# Output files:
#   - Figure_2E_volcano_CTRL_EGS.png
#   - Figure_2E_volcano_CTRL_EGS.pdf
#   - Figure_2E_volcano_CTRL_EGS.svg
#   - Figure_2E_volcano_CTRL_EGS_results_table.csv
#   - Figure_2E_volcano_CTRL_EGS_labelled_proteins.csv
#   - Figure_2E_volcano_CTRL_EGS_sessionInfo.txt
# ============================================================


# ----------------------------
# 0. Load packages
# ----------------------------

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(limma)
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

cat("\nSample mapping used for Figure 2E:\n")
print(mapping_tbl)

colnames(Filtered_final)[match(sample_cols, colnames(Filtered_final))] <- mapping_tbl$new_name


# ----------------------------
# 4. Define final retained sample groups
# ----------------------------

ctrl_keep <- c("CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5")
egs_keep  <- c("EGS1", "EGS2", "EGS3")
keep_samples <- c(ctrl_keep, egs_keep)

stopifnot(all(keep_samples %in% colnames(Filtered_final)))

cat("\nRetained samples used in volcano analysis:\n")
print(keep_samples)

cat("\nDetected proteins per retained sample in Filtered_final.tsv:\n")
print(
  colSums(
    !is.na(Filtered_final[, keep_samples]) &
      is.finite(as.matrix(Filtered_final[, keep_samples])) &
      as.matrix(Filtered_final[, keep_samples]) > 0
  )
)


# ----------------------------
# 5. Build annotation table
# ----------------------------

annot_tbl <- Filtered_final %>%
  transmute(
    Protein = PG.ProteinGroups,
    Gene = PG.Genes,
    Protein_description = PG.ProteinDescriptions
  ) %>%
  mutate(
    Gene = sub(";.*", "", Gene),
    Gene = str_trim(Gene),
    Gene = if_else(is.na(Gene) | Gene == "", Protein, Gene)
  ) %>%
  distinct(Protein, .keep_all = TRUE)


# ----------------------------
# 6. Build protein-level log2 matrix
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
    Protein = PG.ProteinGroups,
    all_of(keep_samples)
  ) %>%
  mutate(
    across(all_of(keep_samples), ~ as.numeric(.))
  ) %>%
  group_by(Protein) %>%
  summarise(
    across(all_of(keep_samples), safe_max),
    .groups = "drop"
  )

mat_raw <- mat_tbl %>%
  column_to_rownames("Protein") %>%
  as.matrix()

mode(mat_raw) <- "numeric"

mat_log2 <- log2(mat_raw + 1)


# ----------------------------
# 7. Require at least 2 valid values per group
# ----------------------------

keep_proteins <- rowSums(is.finite(mat_log2[, ctrl_keep, drop = FALSE])) >= 2 &
  rowSums(is.finite(mat_log2[, egs_keep, drop = FALSE])) >= 2

x <- mat_log2[keep_proteins, keep_samples, drop = FALSE]

cat("\nProteins retained after >=2 valid values per group:", nrow(x), "\n")


# ----------------------------
# 8. Remove zero-variance proteins
# ----------------------------

x <- x[
  apply(x, 1, function(v) sd(v, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

cat("\nProteins entering limma volcano analysis:", nrow(x), "\n")


# ----------------------------
# 9. Row-median imputation for limma
# ----------------------------

x_imp <- t(apply(x, 1, function(v) {
  v[!is.finite(v) | is.na(v)] <- median(v[is.finite(v)], na.rm = TRUE)
  v
}))

colnames(x_imp) <- colnames(x)
rownames(x_imp) <- rownames(x)


# ----------------------------
# 10. Median-normalise per sample
# ----------------------------

sample_medians <- apply(
  x_imp,
  2,
  function(v) median(v[is.finite(v)], na.rm = TRUE)
)

x_norm <- sweep(
  x_imp,
  2,
  sample_medians,
  FUN = "-"
)


# ----------------------------
# 11. limma differential abundance analysis
# ----------------------------

sample_info <- tibble(
  Sample = keep_samples,
  Group = c(rep("CTRL", length(ctrl_keep)), rep("EGS", length(egs_keep)))
) %>%
  mutate(
    Group = factor(Group, levels = c("CTRL", "EGS"))
  )

x_norm <- x_norm[, sample_info$Sample, drop = FALSE]

design <- model.matrix(~ 0 + Group, data = sample_info)
colnames(design) <- levels(sample_info$Group)

fit <- lmFit(x_norm, design)

contrast_matrix <- makeContrasts(
  EGS_vs_CTRL = EGS - CTRL,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

tt <- topTable(
  fit2,
  coef = "EGS_vs_CTRL",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
) %>%
  rownames_to_column("Protein")


# ----------------------------
# 12. Build volcano results table
# ----------------------------

fc_thr <- 0.8
p_thr  <- 0.1

res_plot <- tt %>%
  transmute(
    Protein,
    log2FC = logFC,
    p_value = P.Value,
    adjusted_p_value = adj.P.Val
  ) %>%
  left_join(annot_tbl, by = "Protein") %>%
  mutate(
    p_value = if_else(is.na(p_value) | p_value <= 0, 1e-300, p_value),
    adjusted_p_value = if_else(
      is.na(adjusted_p_value) | adjusted_p_value <= 0,
      1e-300,
      adjusted_p_value
    ),
    neglog10_p = -log10(p_value),
    Category = case_when(
      p_value < p_thr & log2FC >= fc_thr ~ "Higher in EGS",
      p_value < p_thr & log2FC <= -fc_thr ~ "Higher in CTRL",
      TRUE ~ "Not significant"
    ),
    Category = factor(
      Category,
      levels = c("Higher in CTRL", "Higher in EGS", "Not significant")
    )
  )

cat("\nVolcano category counts:\n")
print(table(res_plot$Category))


# ----------------------------
# 13. Select proteins to label
# ----------------------------

n_labels_each_side <- 6

label_egs <- res_plot %>%
  filter(Category == "Higher in EGS") %>%
  arrange(p_value, desc(log2FC)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  slice_head(n = n_labels_each_side)

label_ctrl <- res_plot %>%
  filter(Category == "Higher in CTRL") %>%
  arrange(p_value, log2FC) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  slice_head(n = n_labels_each_side)

label_tbl <- bind_rows(label_ctrl, label_egs) %>%
  distinct(Gene, .keep_all = TRUE)

cat("\nLabelled proteins:\n")
print(
  label_tbl %>%
    select(Gene, log2FC, p_value, adjusted_p_value, Category)
)


# ----------------------------
# 14. Colours
# ----------------------------

pal <- c(
  "Higher in CTRL" = "#15a2a2",
  "Higher in EGS" = "#fea802",
  "Not significant" = "grey82"
)


# ----------------------------
# 15. Axis limits
# ----------------------------

xmax <- ceiling(max(abs(res_plot$log2FC), na.rm = TRUE) * 10) / 10
ymax <- ceiling(max(res_plot$neglog10_p, na.rm = TRUE) * 10) / 10

x_pad <- 0.80
y_pad <- 0.35


# ----------------------------
# 16. Volcano plot
# ----------------------------

p_volcano <- ggplot() +
  geom_point(
    data = filter(res_plot, Category == "Not significant"),
    aes(x = log2FC, y = neglog10_p, colour = Category),
    size = 1.8,
    alpha = 0.55
  ) +
  geom_point(
    data = filter(res_plot, Category == "Higher in CTRL"),
    aes(x = log2FC, y = neglog10_p, colour = Category),
    size = 2.8,
    alpha = 0.95
  ) +
  geom_point(
    data = filter(res_plot, Category == "Higher in EGS"),
    aes(x = log2FC, y = neglog10_p, colour = Category),
    size = 2.8,
    alpha = 0.95
  ) +
  geom_vline(
    xintercept = c(-fc_thr, fc_thr),
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey60"
  ) +
  geom_hline(
    yintercept = -log10(p_thr),
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey60"
  ) +
  geom_label_repel(
    data = filter(label_tbl, log2FC < 0),
    aes(x = log2FC, y = neglog10_p, label = Gene),
    colour = pal["Higher in CTRL"],
    fill = alpha("white", 0.85),
    label.size = NA,
    size = 4.3,
    fontface = "bold",
    direction = "both",
    hjust = 1,
    nudge_x = -0.30,
    box.padding = 0.55,
    point.padding = 0.45,
    segment.color = "grey45",
    segment.size = 0.35,
    min.segment.length = 0,
    max.overlaps = Inf,
    force = 8,
    force_pull = 0.3,
    seed = 123
  ) +
  geom_label_repel(
    data = filter(label_tbl, log2FC > 0),
    aes(x = log2FC, y = neglog10_p, label = Gene),
    colour = pal["Higher in EGS"],
    fill = alpha("white", 0.85),
    label.size = NA,
    size = 4.3,
    fontface = "bold",
    direction = "both",
    hjust = 0,
    nudge_x = 0.30,
    box.padding = 0.55,
    point.padding = 0.45,
    segment.color = "grey45",
    segment.size = 0.35,
    min.segment.length = 0,
    max.overlaps = Inf,
    force = 8,
    force_pull = 0.3,
    seed = 123
  ) +
  scale_colour_manual(
    values = pal,
    drop = FALSE
  ) +
  coord_cartesian(
    xlim = c(-xmax - x_pad, xmax + x_pad),
    ylim = c(0, ymax + y_pad),
    clip = "off"
  ) +
  labs(
    x = expression(log[2] * " fold change (EGS vs control)"),
    y = expression(-log[10] * italic(P)),
    colour = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.ticks = element_line(linewidth = 0.7, colour = "black"),
    plot.margin = margin(8, 35, 8, 8)
  )

print(p_volcano)


# ----------------------------
# 17. Save outputs
# ----------------------------

write_csv(
  res_plot,
  "Figure_2E_volcano_CTRL_EGS_results_table.csv"
)

write_csv(
  label_tbl,
  "Figure_2E_volcano_CTRL_EGS_labelled_proteins.csv"
)

ggsave(
  filename = "Figure_2E_volcano_CTRL_EGS.png",
  plot = p_volcano,
  width = 9.4,
  height = 6.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Figure_2E_volcano_CTRL_EGS.pdf",
  plot = p_volcano,
  width = 9.4,
  height = 6.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "Figure_2E_volcano_CTRL_EGS.svg",
  plot = p_volcano,
  width = 9.4,
  height = 6.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

sink("Figure_2E_volcano_CTRL_EGS_sessionInfo.txt")
sessionInfo()
sink()


cat("\nSaved:\n")
cat("- Figure_2E_volcano_CTRL_EGS.png\n")
cat("- Figure_2E_volcano_CTRL_EGS.pdf\n")
cat("- Figure_2E_volcano_CTRL_EGS.svg\n")
cat("- Figure_2E_volcano_CTRL_EGS_results_table.csv\n")
cat("- Figure_2E_volcano_CTRL_EGS_labelled_proteins.csv\n")
cat("- Figure_2E_volcano_CTRL_EGS_sessionInfo.txt\n")
