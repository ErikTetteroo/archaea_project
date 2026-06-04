local({

#load project
source("R/load_project.R")

#load in relevant data
coverage_table <- read_csv("Data/cleaned_data/coverage_table.csv")
trna <- read_csv("Data/cleaned_data/trna_summary_c.csv")
missing_common_combos <- read_csv("Data/manuscript/missing_common_trnas.csv")
present_rare_combos <- read_csv("Data/manuscript/rare_trnas.csv")


#seperate the ids that contain the standard set
remaining_ids <- setdiff(
  unique(trna$organism_id),
  union(
    unique(present_rare_combos$organism_id),
    unique(missing_common_combos$organism_id)
  )
)

length(unique(trna$organism_id))  #all organisms
length(remaining_ids)             #organisms with the standard set

#check coverage of standard set
standard_set_ids <- trna %>%
  filter(organism_id %in% remaining_ids)

#Filter trna pseudo genes & low score genes unrecognized anticodons
standard_sets_filtered <- standard_set_ids %>%
  filter(
    inf_score > 30,
    is.na(note) | !str_detect(note, "pseudo"),
    !str_detect(anticodon, "N")   
  )

# pick one organism with standard set
standard_id <- standard_sets_filtered %>%
  filter(
    organism_id == unique(standard_sets_filtered$organism_id)[1]
  )

# merge with coverage table
standard_sets_coverage <- standard_id %>%
  inner_join(
    coverage_table,
    by = "anticodon"
  )

# calculate pairing coverage
pairing_counts <- standard_sets_coverage %>%
  
  # remove initiator methionine tRNAs
  filter(tRNA_type != "iMet") %>%
  
  # handle Ile2 modification
  mutate(
    
    pairing = case_when(
      
      # Ile2 reads AUA after modification
      tRNA_type == "Ile2" &
        anticodon == "CAT" &
        codon == "AUA"
      ~ "M2",
      
      # prevent canonical AUG decoding
      tRNA_type == "Ile2" &
        anticodon == "CAT" &
        codon == "AUG"
      ~ NA_character_,
      
      TRUE ~ pairing
    )
  ) %>%
  
  # keep only real pairings
  filter(!is.na(pairing)) %>%
  
  mutate(
    
    pair_type = case_when(
      pairing == "M"  ~ "M",
      pairing == "M2" ~ "M2",
      pairing == "GU" ~ "GUw",
      pairing == "SU" ~ "SUw"
    )
  ) %>%
  
  group_by(
    organism_id,
    codon,
    pair_type
  ) %>%
  
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  
  pivot_wider(
    names_from = pair_type,
    values_from = count,
    values_fill = 0
  )

pairing_counts <- pairing_counts %>%
  mutate(
    CM = case_when(
      M   > 0 ~ "M",
      M2  > 0 ~ "M2",
      GUw > 0 ~ "GU",
      SUw > 0 ~ "SU",
      TRUE ~ NA_character_
    )
  )

standard_set_coverage_table <- pairing_counts %>%
  mutate(
    CM = factor(
      CM,
      levels = c(NA, "SU", "GU","M2", "M")
    )
  )

source("R/plot_functions/coverage_plot_alternative.R")

codon_table <- give_codon_table()

p <- create_alternative_coverage_plot(standard_set_coverage_table,
                                      codon_table,
                                      coverage_table,
                                      title = "codon coverage of the standard set"
)

# save plot
ggsave(
  filename = "Plots/coverage_plot_standard_set.png",
  plot = p,
  width = 12,
  height = 10,
  dpi = 300
)

standard_set_coverage_table <- standard_set_coverage_table[,-1]
write_csv(standard_set_coverage_table,"Data/manuscript/standard_set_coverage_table.csv")

})