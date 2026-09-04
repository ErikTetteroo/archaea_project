local({

source("R/load_project.R")

#load in data
trna <- read_csv("Data/cleaned_data/trna_summary_c.csv")

# Final filtered dataset
trna_filtered <- trna %>%
  filter(
    inf_score > 30,
    is.na(note) | !str_detect(note, "pseudo"),
    !str_detect(anticodon, "N")
  )

# Everything filtered out
trna_low_quality <- trna %>%
  filter(
    inf_score <= 30 |
      (!is.na(note) & str_detect(note, "pseudo")) |
      str_detect(anticodon, "N")
  )

trna_clean <- trna_filtered[,c(6,2,3)]

#How many unique trnas does each organism contain
unique_trna_per_org <- trna_clean %>%
  distinct(organism_id, tRNA_type, anticodon) %>%
  group_by(organism_id) %>%
  summarise(unique_trna = n(), .groups = "drop")

#count how often each trna type is present in each organism
trna_combo_counts <- trna_clean %>%
  distinct(organism_id, tRNA_type, anticodon) %>%  # remove within-organism duplicates
  group_by(tRNA_type, anticodon) %>%
  summarise(n_organisms = n(), .groups = "drop") %>%
  arrange(desc(n_organisms))


#find out which organisms have rare combos 
trna_presence <- trna_clean %>%
  distinct(organism_id, tRNA_type, anticodon)

rare_combos <- trna_combo_counts %>%
  filter(n_organisms <= 21)

rare_combo_organisms <- rare_combos %>%
  inner_join(trna_presence, by = c("tRNA_type", "anticodon")) %>%
  arrange(tRNA_type, anticodon)

# or are missing common combos
common_combos <- trna_combo_counts %>%
  filter(n_organisms >= 500)

all_organisms <- trna_clean %>%
  distinct(organism_id)

missing_common_combos <- common_combos %>%
  select(tRNA_type, anticodon) %>%
  crossing(all_organisms) %>%   
  anti_join(trna_presence,
            by = c("organism_id", "tRNA_type", "anticodon")) %>%
  arrange(tRNA_type, anticodon)

#condense results
rare_summary <- rare_combo_organisms %>%
  group_by(tRNA_type, anticodon) %>%
  summarise(organisms = list(organism_id), .groups = "drop")

missing_summary <- missing_common_combos %>%
  group_by(tRNA_type, anticodon) %>%
  summarise(missing_organisms = list(organism_id), .groups = "drop")

#save a copy of both
write_csv(missing_common_combos,"Data/manuscript/missing_common_trnas.csv")
write_csv(rare_combo_organisms,"Data/manuscript/rare_trnas.csv")

})