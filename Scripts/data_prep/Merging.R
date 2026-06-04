local({
  
#load packages
source("R/load_project.R")

#read in codon usage, trna data, coverage table, and taxonomy
coverage_table <- read_csv("Data/cleaned_data/coverage_table.csv")
codon_usage <- read_csv("Data/Raw_data/combined_codon_usage.csv")
trna_dat <- read_csv("Data/cleaned_data/trna_summary_c.csv")
taxa <- read_tsv("Data/Raw_data/lineages.tsv")

#Filter trna pseudo genes & low score genes unrecognized anticodons
trna_filtered <- trna_dat %>%
  filter(
    inf_score > 30,
    is.na(note) | !str_detect(note, "pseudo"),
    !str_detect(anticodon, "N")   
  )

# Merge trna & codon usage by organism

codon_trna_pairs <- codon_usage %>%
  rename(organism_id = organism) %>%
  inner_join(trna_filtered, by = "organism_id")

# Ignore initiator Met tRNAs

codon_trna_pairs <- codon_trna_pairs %>%
  filter(tRNA_type != "iMet")

# Merge with coverage table

codon_trna_pairs <- codon_trna_pairs %>%
  inner_join(
    coverage_table,
    by = c("codon", "anticodon")
  )

# Special handling for Ile2(CAU)

codon_trna_pairs <- codon_trna_pairs %>%
  mutate(
    
    pairing = case_when(
      
      # Ile2 specifically reads AUA after modification
      tRNA_type == "Ile2" &
        anticodon == "CAT" &
        codon == "AUA"
      ~ "M2",
      
      # prevent canonical AUG decoding by Ile2
      tRNA_type == "Ile2" &
        anticodon == "CAT" &
        codon == "AUG"
      ~ NA_character_,
      
      TRUE ~ pairing
    )
  )

pairing_counts <- codon_trna_pairs %>%
  filter(!is.na(pairing)) %>%
  mutate(
    pair_type = case_when(
      pairing == "M"  ~ "M",
      pairing == "M2" ~ "M2",
      pairing == "GU" ~ "GUw",
      pairing == "SU" ~ "SUw"
    )
  ) %>%
  group_by(organism_id, codon, pair_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = pair_type,
    values_from = count,
    values_fill = 0
  )

final_data <- codon_usage %>%
  rename(organism_id = organism) %>%
  left_join(pairing_counts,
            by = c("organism_id", "codon")) %>%
  mutate(
    M   = coalesce(M, 0),
    M2  = coalesce(M2, 0),
    GUw = coalesce(GUw, 0),
    SUw = coalesce(SUw, 0)
  )

final_data <- final_data %>%
  mutate(
    CM = case_when(
      M   > 0 ~ "M",
      M2  > 0 ~ "M2",
      GUw > 0 ~ "GU",
      SUw > 0 ~ "SU",
      TRUE ~ NA_character_
    )
  )

final_data <- final_data %>%
  mutate(
    CM = factor(
      CM,
      levels = c(NA, "SU", "GU","M2", "M")
    )
  )

final_data <- final_data %>%
  filter(
    !(amino_acid == "TER" & is.na(CM))
  )

#add taxonomy
colnames(taxa)[1] <- "organism_id"
taxac <- taxa[,-3:-4]

final_data_t <- left_join(final_data, taxac, by = "organism_id")

#export merged file
write_csv(final_data_t, 
          "Data/cleaned_data/merged_codon_usage_data.csv")
})