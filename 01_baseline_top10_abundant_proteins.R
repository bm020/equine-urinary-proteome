# ============================================================
# Table 1: Top 10 most abundant proteins in baseline equine urine
# ============================================================
#
# Manuscript:
# Defining the Equine Urinary Proteome: A Reference Baseline
# for Biomarker Discovery
#
# Purpose:
# This script generates the Top 10 most abundant proteins table
# from the baseline fresh equine urine samples U1, U2 and U3.
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
# Processing:
#   - Intensities are log2-transformed using log2(intensity + 1).
#   - Duplicate gene entries are collapsed within each sample by
#     retaining the maximum log2 intensity value.
#   - Proteins are ranked by median log2 intensity across detected samples.
#
# Display label fixes:
#   - LOC100067869 is displayed as HP.
#   - LOC100062179 is displayed as Pepsin A.
#   - Q95182 is retained as Q95182.
#
# Notes:
#   - Raw identifiers are retained in the mapping output.
#   - Display labels are used only for manuscript readability.
#
# Package/method references:
#   - readr::read_tsv for importing TSV files:
#     https://readr.tidyverse.org/reference/read_delim.html
#   - dplyr::group_by and dplyr::summarise for grouped summaries:
#     https://dplyr.tidyverse.org/reference/group_by.html
#     https://dplyr.tidyverse.org/reference/summarise.html
# ============================================================


# ----------------------------
# 0. Load packages
# ----------------------------

library(readr)
library(dplyr)
library(stringr)


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
      Protein.Group,
      Genes,
      Protein_name = `First.Protein.Description`,
      proteotypic_peptides = `N.Proteotypic.Sequences`
    ),
  U2 %>%
    select(
      Protein.Group,
      Genes,
      Protein_name = `First.Protein.Description`,
      proteotypic_peptides = `N.Proteotypic.Sequences`
    )
) %>%
  group_by(Protein.Group) %>%
  summarise(
    Genes = {
      x <- Genes[!is.na(Genes) & Genes != ""]
      if (length(x) == 0) NA_character_ else x[1]
    },
    Protein_name = {
      x <- Protein_name[!is.na(Protein_name) & Protein_name != ""]
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
    ),
    Gene_symbol = sub(";.*", "", Genes),
    Gene_symbol = str_trim(Gene_symbol),
    Gene_symbol = ifelse(
      is.na(Gene_symbol) | Gene_symbol == "",
      Protein.Group,
      Gene_symbol
    ),
    Display_label = case_when(
      Gene_symbol == "LOC100067869" ~ "HP",
      Gene_symbol == "LOC100062179" ~ "Pepsin A",
      Protein.Group == "Q95182" ~ "Q95182",
      Gene_symbol == "Q95182" ~ "Q95182",
      str_detect(
        Protein_name,
        regex("Major allergen Equ c 1|Equ c 1", ignore_case = TRUE)
      ) ~ "Q95182",
      TRUE ~ Gene_symbol
    )
  )


# ----------------------------
# 4. Create sample-specific tables
# ----------------------------

U1_tbl <- U1U3 %>%
  select(Protein.Group, intensity = all_of(u1_col)) %>%
  mutate(Sample = "U1")

U2_tbl <- U2 %>%
  select(Protein.Group, intensity = all_of(u2_col)) %>%
  mutate(Sample = "U2")

U3_tbl <- U1U3 %>%
  select(Protein.Group, intensity = all_of(u3_col)) %>%
  mutate(Sample = "U3")


# ----------------------------
# 5. Combine samples and calculate log2 intensity
# ----------------------------

all_tbl <- bind_rows(U1_tbl, U2_tbl, U3_tbl) %>%
  left_join(annot_tbl, by = "Protein.Group") %>%
  mutate(
    intensity = as.numeric(intensity),
    log2_intensity = log2(intensity + 1)
  ) %>%
  filter(!is.na(log2_intensity), is.finite(log2_intensity))


# ----------------------------
# 6. Collapse duplicate gene entries within each sample
# ----------------------------

gene_tbl <- all_tbl %>%
  filter(!is.na(Gene_symbol), Gene_symbol != "") %>%
  group_by(Sample, Gene_symbol) %>%
  summarise(
    Display_label = first(Display_label),
    Protein_name = first(Protein_name),
    Protein.Group = first(Protein.Group),
    proteotypic_peptides = max(proteotypic_peptides, na.rm = TRUE),
    log2_intensity = max(log2_intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    proteotypic_peptides = ifelse(
      is.infinite(proteotypic_peptides),
      NA_real_,
      proteotypic_peptides
    )
  ) %>%
  group_by(Gene_symbol) %>%
  filter(n_distinct(Sample) >= 2) %>%
  ungroup() %>%
  filter(!is.na(proteotypic_peptides), proteotypic_peptides >= 2)


# ----------------------------
# 7. Rank proteins by median abundance across samples
# ----------------------------

rank_tbl <- gene_tbl %>%
  group_by(Gene_symbol) %>%
  summarise(
    Display_label = first(Display_label),
    Protein_name = first(Protein_name),
    Protein.Group = first(Protein.Group),
    proteotypic_peptides = max(proteotypic_peptides, na.rm = TRUE),
    detected_samples = n_distinct(Sample),
    median_log2_intensity = median(log2_intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_log2_intensity)) %>%
  mutate(rank = row_number())


# ----------------------------
# 8. Extract Top 10 proteins for manuscript table
# ----------------------------

top10_table <- rank_tbl %>%
  slice_head(n = 10) %>%
  transmute(
    `Protein rank` = rank,
    `Gene symbol` = Display_label,
    `Protein name` = Protein_name,
    `Proteotypic peptides` = proteotypic_peptides,
    `Median log2 intensity` = round(median_log2_intensity, 1)
  )


# ----------------------------
# 9. Save full mapping table for traceability
# ----------------------------

top10_mapping <- rank_tbl %>%
  slice_head(n = 10) %>%
  transmute(
    `Protein rank` = rank,
    `Display label` = Display_label,
    `Gene / accession` = Gene_symbol,
    `Protein group` = Protein.Group,
    `Protein name` = Protein_name,
    `Detected samples` = detected_samples,
    `Proteotypic peptides` = proteotypic_peptides,
    `Median log2 intensity` = round(median_log2_intensity, 1)
  )


# ----------------------------
# 10. Display and save
# ----------------------------

print(top10_table, n = 10, width = Inf)

write_csv(
  top10_table,
  "Table1_top10_most_abundant_baseline_equine_urine_proteins.csv"
)

write_csv(
  top10_mapping,
  "Table1_top10_most_abundant_baseline_equine_urine_proteins_mapping.csv"
)

cat("\nSaved:\n")
cat("- Table1_top10_most_abundant_baseline_equine_urine_proteins.csv\n")
cat("- Table1_top10_most_abundant_baseline_equine_urine_proteins_mapping.csv\n")
