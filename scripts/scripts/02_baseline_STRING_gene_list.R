# ============================================================
# Figure 1D-E: Baseline urinary proteome STRING input list
# ============================================================
#
# Manuscript:
# Defining the Equine Urinary Proteome: A Reference Baseline
# for Biomarker Discovery
#
# Purpose:
# This script generates the gene-symbol input list used for
# STRING v12.0 functional enrichment analysis of the baseline
# fresh equine urinary proteome.
#
# Input files:
#   - Horse_urine_U1&U3_horse_fasta_020625_report.pg_matrix.tsv
#   - Horse_urine_U2_donkey_fasta_160525_report.pg_matrix.tsv
#
# Filtering:
#   1. Protein/gene entries are retained if detected in at least
#      2 of 3 fresh urine samples.
#   2. Protein/gene entries are retained if supported by at least
#      2 proteotypic peptides.
#
# STRING settings applied manually on the STRING website:
#   - STRING version: 12.0
#   - Organism: Equus caballus
#   - Background: whole genome
#   - Enrichment categories: Gene Ontology Biological Process
#     and Reactome pathways
#   - FDR threshold: <= 0.05
#   - Minimum gene count: 2
#   - Functional similarity grouping: >= 0.8
#
# Output files:
#   - Baseline_STRING_input_filtered_proteins.csv
#   - Baseline_STRING_gene_symbols.txt
#   - Baseline_STRING_gene_symbols.csv
#
# Package/method references:
#   - STRING: https://string-db.org/
#   - readr::read_tsv: https://readr.tidyverse.org/reference/read_delim.html
#   - dplyr::group_by and dplyr::summarise:
#     https://dplyr.tidyverse.org/reference/group_by.html
#     https://dplyr.tidyverse.org/reference/summarise.html
# ============================================================


# ----------------------------
# 0. Load packages
# ----------------------------

library(readr)
library(dplyr)
library(stringr)
library(tibble)


# ----------------------------
# 1. Read input files
# ----------------------------

u1u3_file <- "Horse_urine_U1&U3_horse_fasta_020625_report.pg_matrix.tsv"
u2_file   <- "Horse_urine_U2_donkey_fasta_160525_report.pg_matrix.tsv"

stopifnot(file.exists(u1u3_file))
stopifnot(file.exists(u2_file))

U1U3 <- read_tsv(u1u3_file, show_col_types = FALSE)
U2   <- read_tsv(u2_file, show_col_types = FALSE)


# ----------------------------
# 2. Identify intensity columns
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
# 3. Build annotation table
# ----------------------------

annot_tbl <- bind_rows(
  U1U3 %>%
    select(
      Protein_ID = Protein.Group,
      Gene = Genes,
      Protein_name = Protein.Names,
      Description = First.Protein.Description,
      proteotypic_peptides = `N.Proteotypic.Sequences`
    ),
  U2 %>%
    select(
      Protein_ID = Protein.Group,
      Gene = Genes,
      Protein_name = Protein.Names,
      Description = First.Protein.Description,
      proteotypic_peptides = `N.Proteotypic.Sequences`
    )
) %>%
  group_by(Protein_ID) %>%
  summarise(
    Gene = {
      x <- Gene[!is.na(Gene) & Gene != ""]
      if (length(x) == 0) NA_character_ else x[1]
    },
    Protein_name = {
      x <- Protein_name[!is.na(Protein_name) & Protein_name != ""]
      if (length(x) == 0) NA_character_ else x[1]
    },
    Description = {
      x <- Description[!is.na(Description) & Description != ""]
      if (length(x) == 0) NA_character_ else x[1]
    },
    proteotypic_peptides = max(as.numeric(proteotypic_peptides), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    proteotypic_peptides = ifelse(
      is.infinite(proteotypic_peptides),
      NA_real_,
      proteotypic_peptides
    )
  )


# ----------------------------
# 4. Create sample-specific intensity tables
# ----------------------------

U1_tbl <- U1U3 %>%
  select(Protein_ID = Protein.Group, intensity = all_of(u1_col)) %>%
  mutate(Sample = "U1")

U2_tbl <- U2 %>%
  select(Protein_ID = Protein.Group, intensity = all_of(u2_col)) %>%
  mutate(Sample = "U2")

U3_tbl <- U1U3 %>%
  select(Protein_ID = Protein.Group, intensity = all_of(u3_col)) %>%
  mutate(Sample = "U3")


# ----------------------------
# 5. Combine samples
# ----------------------------

all_tbl <- bind_rows(U1_tbl, U2_tbl, U3_tbl) %>%
  left_join(annot_tbl, by = "Protein_ID") %>%
  mutate(
    intensity = as.numeric(intensity),
    detected = !is.na(intensity) & is.finite(intensity) & intensity > 0
  )


# ----------------------------
# 6. Apply baseline filtering
# ----------------------------

baseline_filtered <- all_tbl %>%
  filter(detected) %>%
  group_by(Protein_ID) %>%
  summarise(
    Gene = first(Gene[!is.na(Gene) & Gene != ""]),
    Protein_name = first(Protein_name[!is.na(Protein_name) & Protein_name != ""]),
    Description = first(Description[!is.na(Description) & Description != ""]),
    detected_samples = n_distinct(Sample),
    proteotypic_peptides = max(proteotypic_peptides, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Gene = ifelse(length(Gene) == 0, NA_character_, Gene),
    Protein_name = ifelse(length(Protein_name) == 0, NA_character_, Protein_name),
    Description = ifelse(length(Description) == 0, NA_character_, Description),
    proteotypic_peptides = ifelse(
      is.infinite(proteotypic_peptides),
      NA_real_,
      proteotypic_peptides
    )
  ) %>%
  filter(
    detected_samples >= 2,
    !is.na(proteotypic_peptides),
    proteotypic_peptides >= 2
  )

cat("\nBaseline proteins retained after filtering:", nrow(baseline_filtered), "\n")


# ----------------------------
# 7. Clean gene symbols for STRING
# ----------------------------

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

baseline_genes <- baseline_filtered %>%
  pull(Gene) %>%
  clean_gene_symbols()

cat("Baseline gene symbols exported for STRING:", length(baseline_genes), "\n")


# ----------------------------
# 8. Save outputs
# ----------------------------

write_csv(
  baseline_filtered,
  "Baseline_STRING_input_filtered_proteins.csv"
)

writeLines(
  baseline_genes,
  "Baseline_STRING_gene_symbols.txt"
)

write_csv(
  tibble(Gene = baseline_genes),
  "Baseline_STRING_gene_symbols.csv"
)


# ----------------------------
# 9. Preview gene list
# ----------------------------

cat("\nFirst 30 genes in baseline STRING list:\n")
print(head(baseline_genes, 30))

cat("\nSaved:\n")
cat("- Baseline_STRING_input_filtered_proteins.csv\n")
cat("- Baseline_STRING_gene_symbols.txt\n")
cat("- Baseline_STRING_gene_symbols.csv\n")
