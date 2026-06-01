library(tidyverse)
library(dplyr)
library(tidyr)

standard_set_coverage_table <- read_csv("Data/standard_set_coverage_table.csv")
missing_common_combos <- read_csv("Data/missing_common_combos.csv")
present_rare_combos <- read_csv("Data/present_rare_combos.csv")
coverage_table <- read_csv("Data/Coverage_table.csv")
coverage_table <- coverage_table[,-1]

missing_df <- missing_common_combos %>%
  mutate(
    tRNA_type = if_else(tRNA_type == "Ile2", "Ile", tRNA_type),
    change_type = "loss"
  )

present_df <- present_rare_combos %>%
  mutate(
    tRNA_type = if_else(tRNA_type == "Ile2", "Ile", tRNA_type),
    change_type = "gain"
  )

trna_changes <- bind_rows(
  missing_df,
  present_df
)

aa_states <- trna_changes %>%
  mutate(
    event = paste0(
      ifelse(change_type == "loss", "-", "+"),
      anticodon
    )
  ) %>%
  group_by(organism_id, tRNA_type) %>%
  summarise(
    aa_signature = paste(sort(event), collapse = ";"),
    .groups = "drop"
  )

aa_state_counts <- aa_states %>%
  count(tRNA_type, aa_signature, sort = TRUE)

#organism_states <- trna_changes %>%
#  mutate(
#    event = paste0(
#      ifelse(change_type == "loss", "-", "+"),
#      tRNA_type, "_",
#      anticodon
#    )
#  ) %>%
#  group_by(organism_id) %>%
#  summarise(
#    state_signature = paste(sort(event), collapse = ";"),
#    .groups = "drop"
#  )
#
#state_counts <- organism_states %>%
#  count(state_signature, sort = TRUE)


coverage_effects <- trna_changes %>%
  left_join(
    coverage_table,
    by = "anticodon"
  ) %>%
  filter(value != "0")

coverage_effects <- coverage_effects %>%
  mutate(
    delta = case_when(
      value == "M"  & change_type == "gain" ~  1,
      value == "M"  & change_type == "loss" ~ -1,
      value == "GU" & change_type == "gain" ~  0.5,
      value == "GU" & change_type == "loss" ~ -0.5
    )
  )

codon_delta_profiles <- coverage_effects %>%
  group_by(organism_id, tRNA_type, codon) %>%
  summarise(
    delta_M = sum(delta[value == "M"], na.rm = TRUE),
    delta_GUw = sum(delta[value == "GU"], na.rm = TRUE),
    delta_CV = sum(delta),
    .groups = "drop"
  )

aa_delta_state_counts <- codon_delta_profiles %>%
  group_by(
    tRNA_type,
    codon,
    delta_M,
    delta_GUw,
    delta_CV
  ) %>%
  summarise(
    n = n_distinct(organism_id),
    organisms = paste(sort(unique(organism_id)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(tRNA_type, desc(n))

aa_delta_states <- codon_delta_profiles %>%
  
  # Create codon-level event strings
  mutate(
    event = paste0(
      codon,
      "(",
      delta_M, ",",
      delta_GUw, ",",
      delta_CV,
      ")"
    )
  ) %>%
  
  # Collapse into amino-acid-specific signatures
  group_by(organism_id, tRNA_type) %>%
  summarise(
    aa_signature = paste(sort(event), collapse = ";"),
    .groups = "drop"
  )

# Count recurring amino-acid decoding states
aa_delta_state_counts <- aa_delta_states %>%
  count(tRNA_type, aa_signature, sort = TRUE)

#write_csv(aa_delta_state_counts,"Data/aa_deviating_state_counts.csv")
