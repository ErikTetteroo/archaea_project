local({
  
# Archaeal tRNA Coverage Plotting Functions
source("R/load_project.R")

standard_set_coverage_table <- read_csv("Data/manuscript/standard_set_coverage_table.csv")
deviating_aa_sets_per_id <- read_csv("Data/manuscript/aa_deviating_state_id_lookup.csv")
merged <- read_csv("Data/cleaned_data/merged_codon_usage_data.csv")

all_codons <- merged %>%
  distinct(codon)

standard_full <- all_codons %>%
  left_join(standard_set_coverage_table, by = "codon") %>%
  mutate(across(c(M, GUw, SUw, M2), ~replace_na(.x, 0))) %>%
  mutate(
    CM = case_when(
      M  > 0 ~ "M",
      M2 > 0 ~ "M2",
      GUw > 0 ~ "GU",
      SUw > 0 ~ "SU",
      TRUE ~ NA_character_
    )
  )


reconstructed <- deviating_aa_sets_per_id %>%
  group_by(organism_id) %>%
  group_modify(~{
    reconstruct_coverage(
      standard_full,
      .x
    )
  }) %>%
  ungroup()

reconstructed <- reconstructed %>%
  mutate(across(M2, ~as.integer(. > 0)))

merged_bin <- merged %>%
  mutate(across(c(M, GUw, SUw, M2), ~as.integer(. > 0)))

comparison <- merged_bin %>%
  left_join(reconstructed, by = c("organism_id", "codon"), suffix = c("_direct", "_recon"))

print("The 2 methods match in the following way")
print(comparison %>%
  summarise(
    M_match   = mean(M_direct == M_recon, na.rm = TRUE),
    GUw_match = mean(GUw_direct == GUw_recon, na.rm = TRUE),
    SUw_match = mean(SUw_direct == SUw_recon, na.rm = TRUE),
    M2_match  = mean(M2_direct == M2_recon, na.rm = TRUE),
    CM_match  = mean(CM_direct == CM_recon, na.rm = TRUE)
  ))

comparison_flagged <- comparison %>%
  mutate(
    M_diff   = M_direct   != M_recon,
    GUw_diff = GUw_direct != GUw_recon,
    SUw_diff = SUw_direct != SUw_recon,
    M2_diff  = M2_direct  != M2_recon,
    CM_diff  = CM_direct  != CM_recon
  )

mismatches <- comparison_flagged %>%
  filter(
    M_diff | GUw_diff | SUw_diff | M2_diff | CM_diff
  ) %>%
  select(
    organism_id,
    codon,
    amino_acid,
    M_direct, M_recon,
    GUw_direct, GUw_recon,
    SUw_direct, SUw_recon,
    M2_direct, M2_recon,
    CM_direct, CM_recon
  )
write_csv(mismatches,"Data/manuscript/mismatches_between_deviation_and_direct_method.csv")
})