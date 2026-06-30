# ============================================================
# Figure 1C: Rank-abundance profile of baseline equine urine
# ============================================================
#
# Manuscript:
# Defining the Equine Urinary Proteome: A Reference Baseline
# for Biomarker Discovery
#
# Purpose:
# This script generates the rank-abundance plot for the baseline
# fresh equine urine samples U1, U2 and U3.
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
#   - Sample-wise median-centred log2 intensity is used for plotting.
#
# Display label fixes:
#   - LOC100067869 is displayed as HP.
#   - LOC100062179 is displayed as Pepsin A.
#   - Q95182 is retained as Q95182.
#
# Notes:
#   - Raw identifiers are retained in Protein.Group and Gene_symbol.
#   - Display_label is used only for figure readability.
#
# Package/method references:
#   - readr::read_tsv for importing TSV files:
#     https://readr.tidyverse.org/reference/read_delim.html
#   - dplyr::group_by and dplyr::summarise for grouped summaries:
#     https://dplyr.tidyverse.org/reference/group_by.html
#     https://dplyr.tidyverse.org/reference/summarise.html
#   - ggplot2 for figure generation:
#     https://ggplot2.tidyverse.org/
#   - ggsave for exporting figures:
#     https://ggplot2.tidyverse.org/reference/ggsave.html
# ============================================================


# ----------------------------
# 0. Load packages
# ----------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(svglite)


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
# 4. Build sample tables
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
# 6. Collapse duplicate gene entries within sample
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
# 7. Rank proteins by median log2 intensity across samples
# ----------------------------

rank_tbl <- gene_tbl %>%
  group_by(Gene_symbol) %>%
  summarise(
    Display_label = first(Display_label),
    Protein_name = first(Protein_name),
    Protein.Group = first(Protein.Group),
    proteotypic_peptides = max(proteotypic_peptides, na.rm = TRUE),
    median_log2_intensity = median(log2_intensity, na.rm = TRUE),
    detected_samples = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  arrange(desc(median_log2_intensity)) %>%
  mutate(rank = row_number())


# ----------------------------
# 8. Join rank and calculate sample-wise median-centred intensity
# ----------------------------

plot_df <- gene_tbl %>%
  left_join(rank_tbl %>% select(Gene_symbol, rank), by = "Gene_symbol") %>%
  group_by(Sample) %>%
  mutate(
    median_centred_log2 = log2_intensity -
      median(log2_intensity, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    Sample = factor(Sample, levels = c("U1", "U2", "U3"))
  )


# ----------------------------
# 9. Anchor points for labels
# ----------------------------

anchor_tbl <- plot_df %>%
  group_by(Gene_symbol, rank) %>%
  summarise(
    point_y = median(median_centred_log2, na.rm = TRUE),
    .groups = "drop"
  )


# ----------------------------
# 10. Top-abundance labels
# ----------------------------

top_targets <- c(
  "UMOD",
  "HP",
  "ALB",
  "AMBP",
  "NPC2",
  "F2",
  "H9GZT5",
  "HPX",
  "Q95182",
  "CFB"
)

top_label_df <- rank_tbl %>%
  filter(Display_label %in% top_targets) %>%
  left_join(anchor_tbl, by = c("Gene_symbol", "rank")) %>%
  filter(!is.na(point_y)) %>%
  mutate(display_label = Display_label) %>%
  arrange(desc(point_y))

top_label_df$label_x <- 190

top_label_df$label_y <- seq(
  from = 10.8,
  to   = 3.7,
  length.out = nrow(top_label_df)
)


# ----------------------------
# 11. Low-abundance labels
# ----------------------------

bottom_label_df <- rank_tbl %>%
  left_join(anchor_tbl, by = c("Gene_symbol", "rank")) %>%
  filter(
    !is.na(Gene_symbol),
    Gene_symbol != "",
    !is.na(point_y),
    detected_samples >= 2,
    proteotypic_peptides >= 2,
    rank >= quantile(rank, 0.72, na.rm = TRUE)
  ) %>%
  arrange(median_log2_intensity, desc(point_y)) %>%
  slice_head(n = 10) %>%
  arrange(desc(point_y)) %>%
  mutate(display_label = Display_label)

bottom_label_df$label_x <- 960

bottom_label_df$label_y <- seq(
  from = 2.1,
  to   = -5.1,
  length.out = nrow(bottom_label_df)
)


# ----------------------------
# 12. Colours
# ----------------------------

palette_vals <- c(
  "U1" = "#4f5bff",
  "U2" = "#e5c400",
  "U3" = "#45d15a"
)


# ----------------------------
# 13. Theme
# ----------------------------

pub_theme <- theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 17, colour = "black"),
    axis.text = element_text(size = 13.5, colour = "black"),
    legend.title = element_text(size = 15, colour = "black"),
    legend.text = element_text(size = 13.5, colour = "black"),
    legend.position = "right",
    legend.justification = "center",
    legend.box = "vertical",
    legend.background = element_blank(),
    axis.line = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks = element_line(linewidth = 0.6, colour = "black"),
    plot.margin = margin(12, 25, 12, 18)
  )


# ----------------------------
# 14. Plot
# ----------------------------

p <- ggplot(plot_df, aes(x = rank, y = median_centred_log2, colour = Sample)) +
  geom_point(size = 2.0, alpha = 0.42, stroke = 0) +
  
  geom_segment(
    data = top_label_df,
    aes(
      x = rank,
      xend = label_x - 8,
      y = point_y,
      yend = label_y
    ),
    inherit.aes = FALSE,
    colour = "grey35",
    linewidth = 0.35,
    lineend = "round"
  ) +
  
  geom_text(
    data = top_label_df,
    aes(
      x = label_x,
      y = label_y,
      label = display_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.7,
    colour = "black"
  ) +
  
  geom_segment(
    data = bottom_label_df,
    aes(
      x = rank,
      xend = label_x - 8,
      y = point_y,
      yend = label_y
    ),
    inherit.aes = FALSE,
    colour = "grey35",
    linewidth = 0.35,
    lineend = "round"
  ) +
  
  geom_text(
    data = bottom_label_df,
    aes(
      x = label_x,
      y = label_y,
      label = display_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.7,
    colour = "black"
  ) +
  
  scale_colour_manual(values = palette_vals) +
  
  scale_x_continuous(
    breaks = c(0, 250, 500, 750, 1000),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_y_continuous(
    breaks = c(-5, 0, 5, 10)
  ) +
  
  coord_cartesian(
    xlim = c(0, 1100),
    ylim = c(-6, 11.8),
    clip = "off"
  ) +
  
  labs(
    x = "Protein rank",
    y = expression("Median-centred log"[2] * " intensity"),
    colour = "Sample"
  ) +
  
  guides(
    colour = guide_legend(
      override.aes = list(size = 5, alpha = 1)
    )
  ) +
  
  pub_theme

print(p)


# ----------------------------
# 15. Print labelled proteins for checking
# ----------------------------

cat("\nTop-abundance labelled proteins:\n")
print(
  top_label_df %>%
    select(
      Gene_symbol,
      Display_label,
      Protein.Group,
      Protein_name,
      rank,
      median_log2_intensity,
      point_y
    )
)

cat("\nLow-abundance labelled proteins:\n")
print(
  bottom_label_df %>%
    select(
      Gene_symbol,
      Display_label,
      Protein.Group,
      Protein_name,
      rank,
      median_log2_intensity,
      point_y
    )
)


# ----------------------------
# 16. Save mapping table for reproducibility
# ----------------------------

label_mapping_table <- bind_rows(
  top_label_df %>%
    mutate(Label_category = "Top abundance"),
  bottom_label_df %>%
    mutate(Label_category = "Low abundance")
) %>%
  select(
    Label_category,
    Gene_symbol,
    Display_label,
    Protein.Group,
    Protein_name,
    proteotypic_peptides,
    detected_samples,
    rank,
    median_log2_intensity
  ) %>%
  distinct()

write_csv(
  label_mapping_table,
  "Figure_1C_rank_abundance_label_mapping.csv"
)


# ----------------------------
# 17. Save figure
# ----------------------------

ggsave(
  "Figure_1C_rank_abundance_baseline_urine_final.svg",
  plot = p,
  width = 10.5,
  height = 6.0,
  units = "in"
)

ggsave(
  "Figure_1C_rank_abundance_baseline_urine_final.pdf",
  plot = p,
  width = 10.5,
  height = 6.0,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  "Figure_1C_rank_abundance_baseline_urine_final.png",
  plot = p,
  width = 10.5,
  height = 6.0,
  units = "in",
  dpi = 600
)

cat("\nSaved:\n")
cat("- Figure_1C_rank_abundance_baseline_urine_final.svg\n")
cat("- Figure_1C_rank_abundance_baseline_urine_final.pdf\n")
cat("- Figure_1C_rank_abundance_baseline_urine_final.png\n")
cat("- Figure_1C_rank_abundance_label_mapping.csv\n")
