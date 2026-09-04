#load project
source("R/load_project.R")

# Load data
codon_pair_usage <- read_csv(
  "Data/Raw_data/combined_codon_pair_usage.csv"
)

merged <- read_csv(
  "Data/cleaned_data/merged_codon_usage_data.csv"
)

# Standardize codon format
merged <- merged %>%
  mutate(codon = rna_to_dna(codon))


# Create lookup tables
gc_lookup <- merged %>%
  select(organism_id, GC3) %>%
  distinct()

codon_lookup <- merged %>%
  select(
    organism_id,
    codon,
    CM,
    RSCU
  )

# Add organism-level GC3
codon_pair_usage <- codon_pair_usage %>%
  left_join(
    gc_lookup,
    by = c("organism" = "organism_id")
  )

# Add coverage and RSCU for codon 1
codon_pair_usage <- codon_pair_usage %>%
  left_join(
    codon_lookup,
    by = c(
      "organism" = "organism_id",
      "c1" = "codon"
    )
  ) %>%
  rename(
    CM1 = CM,
    RSCU1 = RSCU
  )

# Add coverage and RSCU for codon 2
codon_pair_usage <- codon_pair_usage %>%
  left_join(
    codon_lookup,
    by = c(
      "organism" = "organism_id",
      "c2" = "codon"
    )
  ) %>%
  rename(
    CM2 = CM,
    RSCU2 = RSCU
  )

# Create pair-level variables
codon_pair_usage <- codon_pair_usage %>%
  mutate(
    CMP = paste0(CM1, "_X_", CM2),
    log2_RDCU = log2(RDCU)
  )

# Save processed dataset
write_csv(
  codon_pair_usage,
  "Data/cleaned_data/codon_pair_usage_clean.csv"
)