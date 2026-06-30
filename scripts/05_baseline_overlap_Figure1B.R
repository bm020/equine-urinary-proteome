# ============================================================
# Figure 1B: Overlap of proteins identified across baseline urine samples
# ============================================================
#
# Manuscript:
# Defining the Equine Urinary Proteome: A Reference Baseline
# for Biomarker Discovery
#
# Purpose:
# This script generates protein lists for the baseline overlap
# analysis of fresh equine urine samples U1, U2 and U3.
#
# Input files:
#   - Horse_urine_U1&U3_horse_fasta_020625_report.pg_matrix.tsv
#   - Horse_urine_U2_donkey_fasta_160525_report.pg_matrix.tsv
#
# Output files:
#   - Figure_1B_U1_protein_groups.txt
#   - Figure_1B_U2_protein_groups.txt
#   - Figure_1B_U3_protein_groups.txt
#   - Figure_1B_baseline_overlap_counts.csv
#
# Notes:
#   - Protein.Group is used as the primary identifier for overlap.
#   - Intensity columns are selected by sample name rather than column order.
#   - Pooled samples are excluded.
# ============================================================


# ----------------------------
# 0. Load packages
# ----------------------------

library(readr)
library(dplyr)
library(stringr)
library(tibble)


# ----------------------------
# 1. Read files
# ----------------------------

u1u3_file <- "Horse_urine_U1&U3_horse_fasta_020625_report.pg_matrix.tsv"
u2_file   <- "Horse_urine_U2_donkey_fasta_160525_report.pg_matrix.tsv"

stopifnot(file.exists(u1u3_file))
stopifnot(file.exists(u2_file))

U1U3 <- read_tsv(u1u3_file, show_col_types = FALSE)
U2   <- read_tsv(u2_file, show_col_types = FALSE)


# ----------------------------
# 2. Identify intensity columns robustly
# ----------------------------

u1u3_candidates <- grep("\\.wiff$", colnames(U1U3), value = TRUE)
u2_candidates   <- grep("\\.wiff$", colnames(U2), value = TRUE)

u1_col <- u1u3_candidates[
  str_detect(u1u3_candidates, "U1_STRAP") &
    !str_detect(u1u3_candidates, regex("Pooled", ignore_case = TRUE))
]

u3_col <- u1u3_candidates[
  str_detect(u1u3_candidates, "U3_STRAP") &
    !str_detect(u1u3_candidates, regex("Pooled", ignore_case = TRUE))
]

u2_col <- u2_candidates[
  str_detect(u2_candidates, "U2_STRAP") &
    !str_detect(u2_candidates, regex("Pooled", ignore_case = TRUE))
]

if (length(u1_col) != 1) stop("Could not uniquely identify the U1 intensity column.")
if (length(u2_col) != 1) stop("Could not uniquely identify the U2 intensity column.")
if (length(u3_col) != 1) stop("Could not uniquely identify the U3 intensity column.")

cat("\nSelected intensity columns:\n")
cat("U1:", u1_col, "\n")
cat("U2:", u2_col, "\n")
cat("U3:", u3_col, "\n")


# ----------------------------
# 3. Build sample-specific detected protein lists
# ----------------------------

get_detected_proteins <- function(df, intensity_col, sample_name) {
  df %>%
    select(
      Protein_ID = Protein.Group,
      Gene = Genes,
      intensity = all_of(intensity_col)
    ) %>%
    mutate(
      Sample = sample_name,
      intensity = as.numeric(intensity),
      detected = !is.na(intensity) & is.finite(intensity) & intensity > 0
    ) %>%
    filter(detected, !is.na(Protein_ID), Protein_ID != "") %>%
    distinct(Sample, Protein_ID, .keep_all = TRUE)
}

U1_tbl <- get_detected_proteins(U1U3, u1_col, "U1")
U2_tbl <- get_detected_proteins(U2,   u2_col, "U2")
U3_tbl <- get_detected_proteins(U1U3, u3_col, "U3")

all_detected <- bind_rows(U1_tbl, U2_tbl, U3_tbl)


# ----------------------------
# 4. Extract unique protein groups per sample
# ----------------------------

proteins_U1 <- all_detected %>%
  filter(Sample == "U1") %>%
  pull(Protein_ID) %>%
  unique()

proteins_U2 <- all_detected %>%
  filter(Sample == "U2") %>%
  pull(Protein_ID) %>%
  unique()

proteins_U3 <- all_detected %>%
  filter(Sample == "U3") %>%
  pull(Protein_ID) %>%
  unique()


# ----------------------------
# 5. Calculate overlap counts
# ----------------------------

overlap_counts <- tibble(
  Category = c(
    "U1",
    "U2",
    "U3",
    "U1_U2",
    "U1_U3",
    "U2_U3",
    "U1_U2_U3"
  ),
  Count = c(
    length(proteins_U1),
    length(proteins_U2),
    length(proteins_U3),
    length(intersect(proteins_U1, proteins_U2)),
    length(intersect(proteins_U1, proteins_U3)),
    length(intersect(proteins_U2, proteins_U3)),
    length(Reduce(intersect, list(proteins_U1, proteins_U2, proteins_U3)))
  )
)

print(overlap_counts)


# ----------------------------
# 6. Save files for Venny or other Venn tools
# ----------------------------

writeLines(
  sort(proteins_U1),
  "Figure_1B_U1_protein_groups.txt"
)

writeLines(
  sort(proteins_U2),
  "Figure_1B_U2_protein_groups.txt"
)

writeLines(
  sort(proteins_U3),
  "Figure_1B_U3_protein_groups.txt"
)

write_csv(
  overlap_counts,
  "Figure_1B_baseline_overlap_counts.csv"
)


# ----------------------------
# 7. Print summary
# ----------------------------

cat("\nProtein groups detected:\n")
cat("U1:", length(proteins_U1), "\n")
cat("U2:", length(proteins_U2), "\n")
cat("U3:", length(proteins_U3), "\n")
cat("Shared across U1, U2 and U3:", length(Reduce(intersect, list(proteins_U1, proteins_U2, proteins_U3))), "\n")

cat("\nSaved:\n")
cat("- Figure_1B_U1_protein_groups.txt\n")
cat("- Figure_1B_U2_protein_groups.txt\n")
cat("- Figure_1B_U3_protein_groups.txt\n")
cat("- Figure_1B_baseline_overlap_counts.csv\n")
