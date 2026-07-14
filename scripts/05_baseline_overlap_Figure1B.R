# ============================================================
# Figure 1B: Baseline urinary proteome overlap by gene symbol
# ============================================================
#
# Purpose:
#   Extract gene symbols detected in each baseline urine sample
#   and export gene lists for Venny/overlap visualisation.
#
# Inputs:
#   - Horse_urine_U1&U3_horse_fasta_020625_report.pg_matrix.tsv
#   - Horse_urine_U2_donkey_fasta_160525_report.pg_matrix.tsv
#
# Outputs:
#   - results/figure1_baseline_overlap/U1_genes.txt
#   - results/figure1_baseline_overlap/U2_genes.txt
#   - results/figure1_baseline_overlap/U3_genes.txt
#   - results/figure1_baseline_overlap/baseline_overlap_summary.csv
#
# Notes:
#   U1 and U3 were analysed together in one Spectronaut/DIA-NN
#   protein-group matrix, while U2 was analysed separately.
#   Gene symbols are used for Figure 1B to support visualisation
#   of shared and sample-specific baseline urinary proteins.
#
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(tibble)

# ============================================================
# 1. File paths
# ============================================================

baseline_u1u3_file <- "data/Horse_urine_U1&U3_horse_fasta_020625_report.pg_matrix.tsv"
baseline_u2_file   <- "data/Horse_urine_U2_donkey_fasta_160525_report.pg_matrix.tsv"

out_dir <- "results/figure1_baseline_overlap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Helper functions
# ============================================================

check_file_exists <- function(file) {
  if (!file.exists(file)) {
    stop(paste("File not found:", file))
  }
}

find_intensity_column <- function(columns, sample_pattern) {
  matched_col <- columns[str_detect(columns, sample_pattern)]
  
  if (length(matched_col) == 0) {
    stop(paste("No intensity column found for pattern:", sample_pattern))
  }
  
  if (length(matched_col) > 1) {
    message("Multiple columns matched pattern: ", sample_pattern)
    message("Using first matched column: ", matched_col[1])
  }
  
  matched_col[1]
}

extract_detected_genes <- function(df, intensity_col) {
  
  required_cols <- c("Protein.Group", "Genes", intensity_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing required column(s):",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  df %>%
    transmute(
      Protein.Group = Protein.Group,
      Genes = Genes,
      intensity = suppressWarnings(as.numeric(.data[[intensity_col]]))
    ) %>%
    filter(!is.na(intensity), intensity > 0) %>%
    mutate(
      Gene = as.character(Genes),
      Gene = sub(";.*", "", Gene),
      Gene = str_trim(Gene),
      Gene = ifelse(Gene == "" | is.na(Gene), NA_character_, Gene)
    ) %>%
    filter(!is.na(Gene)) %>%
    distinct(Gene) %>%
    arrange(Gene) %>%
    pull(Gene)
}

# ============================================================
# 3. Read input files
# ============================================================

check_file_exists(baseline_u1u3_file)
check_file_exists(baseline_u2_file)

baseline_u1u3 <- read_tsv(baseline_u1u3_file, show_col_types = FALSE)
baseline_u2   <- read_tsv(baseline_u2_file, show_col_types = FALSE)

# ============================================================
# 4. Identify sample intensity columns
# ============================================================

u1u3_intensity_cols <- grep("\\.wiff$", colnames(baseline_u1u3), value = TRUE)
u2_intensity_cols   <- grep("\\.wiff$", colnames(baseline_u2), value = TRUE)

u1_col <- find_intensity_column(u1u3_intensity_cols, "U1_STRAP")
u3_col <- find_intensity_column(u1u3_intensity_cols, "U3_STRAP")
u2_col <- find_intensity_column(u2_intensity_cols, "U2_STRAP")

message("Intensity column used for U1: ", u1_col)
message("Intensity column used for U2: ", u2_col)
message("Intensity column used for U3: ", u3_col)

# ============================================================
# 5. Extract detected gene symbols per sample
# ============================================================

genes_u1 <- extract_detected_genes(baseline_u1u3, u1_col)
genes_u2 <- extract_detected_genes(baseline_u2,   u2_col)
genes_u3 <- extract_detected_genes(baseline_u1u3, u3_col)

# ============================================================
# 6. Calculate overlap statistics
# ============================================================

shared_all <- Reduce(intersect, list(genes_u1, genes_u2, genes_u3))
union_all  <- unique(c(genes_u1, genes_u2, genes_u3))

u1_only <- setdiff(genes_u1, union(genes_u2, genes_u3))
u2_only <- setdiff(genes_u2, union(genes_u1, genes_u3))
u3_only <- setdiff(genes_u3, union(genes_u1, genes_u2))

u1_u2_only <- setdiff(intersect(genes_u1, genes_u2), genes_u3)
u1_u3_only <- setdiff(intersect(genes_u1, genes_u3), genes_u2)
u2_u3_only <- setdiff(intersect(genes_u2, genes_u3), genes_u1)

overlap_summary <- tibble(
  Region = c(
    "U1",
    "U2",
    "U3",
    "Union",
    "Shared across U1, U2 and U3",
    "U1 only",
    "U2 only",
    "U3 only",
    "U1 and U2 only",
    "U1 and U3 only",
    "U2 and U3 only"
  ),
  Gene_count = c(
    length(genes_u1),
    length(genes_u2),
    length(genes_u3),
    length(union_all),
    length(shared_all),
    length(u1_only),
    length(u2_only),
    length(u3_only),
    length(u1_u2_only),
    length(u1_u3_only),
    length(u2_u3_only)
  )
)

# ============================================================
# 7. Print summary
# ============================================================

cat("\n===== Figure 1B baseline overlap summary =====\n")
print(overlap_summary)

# ============================================================
# 8. Export gene lists for Venny
# ============================================================

writeLines(sort(unique(genes_u1)), file.path(out_dir, "U1_genes.txt"))
writeLines(sort(unique(genes_u2)), file.path(out_dir, "U2_genes.txt"))
writeLines(sort(unique(genes_u3)), file.path(out_dir, "U3_genes.txt"))

writeLines(sort(unique(shared_all)), file.path(out_dir, "shared_U1_U2_U3_genes.txt"))
writeLines(sort(unique(union_all)),  file.path(out_dir, "union_U1_U2_U3_genes.txt"))

write_csv(
  overlap_summary,
  file.path(out_dir, "baseline_overlap_summary.csv")
)

# ============================================================
# 9. Session information
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo_Figure1B_baseline_overlap.txt")
)
