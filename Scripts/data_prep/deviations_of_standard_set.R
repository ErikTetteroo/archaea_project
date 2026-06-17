local({
  #load project
source("R/load_project.R")

standard_set_coverage_table <- read_csv("Data/manuscript/standard_set_coverage_table.csv")
missing_common_combos <- read_csv("Data/manuscript/missing_common_trnas.csv")
present_rare_combos <- read_csv("Data/manuscript/rare_trnas.csv")
coverage_table <- read_csv("Data/cleaned_data/coverage_table.csv")

#combine the missing common and rare dataframes to form a 'deviations' dataframe
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

#remove iMet changes
trna_changes <- trna_changes %>%
  filter(tRNA_type != "iMet")


#see how many deviating states there are at the organism level
organism_states <- trna_changes %>%
  mutate(
    event = paste0(
      ifelse(change_type == "loss", "-", "+"),
     tRNA_type, "_",
      anticodon
    )
  ) %>%
  group_by(organism_id) %>%
  summarise(
    state_signature = paste(sort(event), collapse = ";"),
    .groups = "drop"
  )

organism_state_counts <- organism_states %>%
  count(state_signature, sort = TRUE)

# check what the effect of deviating states is on codon coverage
coverage_effects <- trna_changes %>%
  left_join(
    coverage_table,
    by = "anticodon"
  ) %>%
  filter(!is.na(pairing))

coverage_effects <- coverage_effects %>%
  mutate(
    delta = case_when(
      pairing == "M"  & change_type == "gain" ~  1,
      pairing == "M"  & change_type == "loss" ~ -1,
      pairing == "GU" & change_type == "gain" ~  1,
      pairing == "GU" & change_type == "loss" ~ -1,
      pairing == "SU" & change_type == "gain" ~  1,
      pairing == "SU" & change_type == "loss" ~ -1,
      pairing == "M2" & change_type == "gain" ~  1,
      pairing == "M2" & change_type == "loss" ~ -1,
      
    )
  )

codon_delta_profiles <- coverage_effects %>%
  group_by(organism_id, tRNA_type, codon) %>%
  summarise(
    delta_M = sum(delta[pairing == "M"], na.rm = TRUE),
    delta_GUw = sum(delta[pairing == "GU"], na.rm = TRUE),
    delta_SUw = sum(delta[pairing == "SU"], na.rm = TRUE),
    delta_M2 = sum(delta[pairing == "M2"], na.rm = TRUE),
    .groups = "drop"
  )

#group by amino acid
aa_delta_states <- codon_delta_profiles %>%
  
  # Create codon-level event strings
  mutate(
    event = paste0(
      codon,
      "(",
      delta_M, ",",
      delta_GUw, ",",
      delta_SUw, ",",
      delta_M2,
      ")"
    )
  ) %>%
  
  # Collapse into amino-acid-specific signatures
  group_by(organism_id, tRNA_type) %>%
  summarise(
    aa_signature = paste(sort(event), collapse = ";"),
    .groups = "drop"
  )

aa_delta_state_counts <- aa_delta_states %>%
  group_by(tRNA_type, aa_signature) %>%
  summarise(
    n = n(),
    
    .groups = "drop"
  ) %>%
  arrange(desc(n))

# output
write_csv(aa_delta_state_counts,"Data/manuscript/aa_deviating_state_counts.csv")
write_csv(aa_delta_states,"Data/manuscript/aa_deviating_state_id_lookup.csv")
write_csv(organism_state_counts,"Data/manuscript/org_deviating_state_counts.csv")
})