# ============================================================
# Figure 3: Functional enrichment input lists for EGS vs CTRL
# ============================================================
#
# Purpose:
# Generate gene-symbol input lists for STRING v12.0 functional
# enrichment analysis of proteins with higher or lower abundance
# in EGS urine compared with CTRL urine.
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
# Selection criteria:
#   - log2FC >= 0.8 and p-value < 0.1: higher in EGS
#   - log2FC <= -0.8 and p-value < 0.1: higher in CTRL
#
# STRING settings applied manually:
#   - STRING version: 12.0
#   - Organism: Equus caballus
#   - Background: whole genome
#   - Enrichment categories: Gene Ontology Biological Process
#     and Reactome pathways
#   - FDR threshold: <= 0.05
#   - Minimum gene count: 2
#   - Functional similarity grouping: >= 0.8
#
# Package/method references:
#   - readr::read_tsv:
#     https://readr.tidyverse.org/reference/read_delim.html
#   - dplyr data manipulation:
#     https://dplyr.tidyverse.org/
#   - stringr string handling:
#     https://stringr.tidyverse.org/
#   - limma empirical Bayes linear modelling:
#     https://bioconductor.org/packages/release/bioc/html/limma.html
#   - STRING functional enrichment:
#     https://string-db.org/
# ============================================================

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(limma)

# ----------------------------
# 1. Read filtered file
# ----------------------------

input_file <- "Filtered_final.tsv"
stopifnot(file.exists(input_file))

Filtered_final <- read_tsv(input_file, show_col_types = FALSE)

# ----------------------------
# 2. Detect and rename sample columns
# ----------------------------

sample_cols <- grep(
  "Horse_urine_U.*PG\\.Quantity$",
  colnames(Filtered_final),
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

cat("\nSample mapping used for Figure 3 enrichment lists:\n")
print(mapping_tbl)

colnames(Filtered_final)[match(sample_cols, colnames(Filtered_final))] <- mapping_tbl$new_name

ctrl_keep <- c("CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5")
egs_keep  <- c("EGS1", "EGS2", "EGS3")
keep_samples <- c(ctrl_keep, egs_keep)

stopifnot(all(keep_samples %in% colnames(Filtered_final)))

# ----------------------------
# 3. Build annotation table
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
# 4. Build log2 protein matrix
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
# 5. Keep proteins with at least 2 valid values per group
# ----------------------------

keep_proteins <- rowSums(is.finite(mat_log2[, ctrl_keep, drop = FALSE])) >= 2 &
  rowSums(is.finite(mat_log2[, egs_keep, drop = FALSE])) >= 2

x <- mat_log2[keep_proteins, keep_samples, drop = FALSE]

x <- x[
  apply(x, 1, function(v) sd(v, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

cat("\nProteins entering limma analysis:", nrow(x), "\n")

# ----------------------------
# 6. Row-median imputation and median normalisation
# ----------------------------

x_imp <- t(apply(x, 1, function(v) {
  v[!is.finite(v) | is.na(v)] <- median(v[is.finite(v)], na.rm = TRUE)
  v
}))

colnames(x_imp) <- colnames(x)
rownames(x_imp) <- rownames(x)

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
# 7. limma differential abundance analysis
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

res <- topTable(
  fit2,
  coef = "EGS_vs_CTRL",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
) %>%
  rownames_to_column("Protein") %>%
  transmute(
    Protein,
    log2FC = logFC,
    p_value = P.Value,
    adjusted_p_value = adj.P.Val
  ) %>%
  left_join(annot_tbl, by = "Protein")

# ----------------------------
# 8. Create STRING input gene lists
# ----------------------------

fc_thr <- 0.8
p_thr  <- 0.1

clean_gene_symbols <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all(";", ",") %>%
    str_split(",") %>%
    unlist() %>%
    str_trim() %>%
    .[!is.na(.)] %>%
    .[. != ""] %>%
    unique()
}

egs_higher_genes <- res %>%
  filter(log2FC >= fc_thr, p_value < p_thr) %>%
  pull(Gene) %>%
  clean_gene_symbols()

ctrl_higher_genes <- res %>%
  filter(log2FC <= -fc_thr, p_value < p_thr) %>%
  pull(Gene) %>%
  clean_gene_symbols()

cat("\nSTRING input list sizes:\n")
cat("Higher abundance in EGS:", length(egs_higher_genes), "\n")
cat("Higher abundance in CTRL:", length(ctrl_higher_genes), "\n")

# ----------------------------
# 9. Save outputs
# ----------------------------

write_csv(
  res,
  "Figure_3_limma_results_for_STRING_enrichment.csv"
)

writeLines(
  egs_higher_genes,
  "Figure_3_STRING_EGS_higher_gene_symbols.txt"
)

writeLines(
  ctrl_higher_genes,
  "Figure_3_STRING_CTRL_higher_gene_symbols.txt"
)

write_csv(
  tibble(Gene = egs_higher_genes),
  "Figure_3_STRING_EGS_higher_gene_symbols.csv"
)

write_csv(
  tibble(Gene = ctrl_higher_genes),
  "Figure_3_STRING_CTRL_higher_gene_symbols.csv"
)

sink("Figure_3_functional_enrichment_input_sessionInfo.txt")
sessionInfo()
sink()

cat("\nSaved:\n")
cat("- Figure_3_limma_results_for_STRING_enrichment.csv\n")
cat("- Figure_3_STRING_EGS_higher_gene_symbols.txt\n")
cat("- Figure_3_STRING_CTRL_higher_gene_symbols.txt\n")
cat("- Figure_3_STRING_EGS_higher_gene_symbols.csv\n")
cat("- Figure_3_STRING_CTRL_higher_gene_symbols.csv\n")
cat("- Figure_3_functional_enrichment_input_sessionInfo.txt\n")
