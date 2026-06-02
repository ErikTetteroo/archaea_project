source("R/load_project.R")

coverage_table <- read_csv("Data/Coverage_table.csv")
trna <- read_csv("Data/Raw_data/trnas_summary.csv")
missing_common_combos <- read_csv("Data/missing_common_combos.csv")
present_rare_combos <- read_csv("Data/present_rare_combos.csv")

#Create consistent organism id
trna <- trna %>%
  mutate(organism_id = str_extract(organism, "GCF_[0-9]+\\.[0-9]+"))

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


sort(unique(missing_common_combos$tRNA_type))==sort(unique(standard_sets_filtered$tRNA_type))
sort(unique(present_rare_combos$tRNA_type))

test <- standard_sets_filtered[standard_sets_filtered$organism_id==unique(standard_sets_filtered$organism_id)[1],]

#Merge with coverage table by codon & anticodon
standard_sets_coverage <- test %>%
  inner_join(coverage_table, by = c("anticodon"))


pairing_counts <- standard_sets_coverage %>%
  
  # Remove initiator methionine tRNAs
  filter(tRNA_type != "iMet") %>%
  
  # Reassign Ile2 AUG pairing to AUA
  mutate(
    codon = case_when(
      tRNA_type == "Ile2" & codon == "AUG" ~ "AUA",
      TRUE ~ codon
    ),
    
    # Optional flag
    ile2_AUA = case_when(
      tRNA_type == "Ile2" & codon == "AUA" ~ TRUE,
      TRUE ~ FALSE
    ),
    
    pair_type = case_when(
      value == "M"  ~ "M",
      value == "GU" ~ "GUw"
    )
  ) %>%
  
  filter(value != "0") %>%
  
  group_by(organism_id, codon, pair_type) %>%
  summarise(
    count = n(),
    ile2_AUA = any(ile2_AUA),
    .groups = "drop"
  ) %>%
  
  pivot_wider(
    names_from = pair_type,
    values_from = count,
    values_fill = 0
  )


pairing_counts$codon
sort(unique(standard_sets_coverage$codon))

#combine pairing counts in Coverage value
GUweight <- 0.5
standard_set_coverage_table <- pairing_counts %>%
  mutate(
    CV = M + GUweight * GUw
  )

standard_set_coverage_table <- standard_set_coverage_table[,-1]
write_csv(standard_set_coverage_table,"Data/standard_set_coverage_table.csv")
