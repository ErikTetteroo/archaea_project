library(tidyverse)
library(dplyr)
library(tidyr)

manuscript_data <- read_csv("Data/Raw_data/archaea_complete_curated_zenodo.csv")
trna_dat <- read.csv("Data/Raw_data/trnas_summary.csv")
matches <- read.csv("Data/manuscript_match_id.csv")

#Create consistent organism id
trna <- trna_dat %>%
  mutate(organism_id = str_extract(organism, "GCF_[0-9]+\\.[0-9]+"))
unique(trna$organism_id)

# Add matching Assembly_Acc to trna
trna_trimmed <- trna %>%
  inner_join(
    matches,
    by = c("organism_id" = "unique.merged.organism_id.")
  )

# Add matching organism_id to manuscript_data
manuscript_trimmed <- manuscript_data %>%
  inner_join(
    matches,
    by = "Assembly_Acc"
  )

manuscript_trimmed_filtered <- manuscript_trimmed %>%
  filter(is.na(curation) | curation == "ok")

which(trna$organism_id=='GCF_000195935.2')

manuscript_trimmed_filtered$tRNA_Type[manuscript_trimmed_filtered$comment_Peter=="splicing went wrong, this is no TTT but CGA and is OK"] <- "Ser"


# -----------------------------
# Clean your trnascan-se results
# -----------------------------
trna_clean <- trna_trimmed %>%
  transmute(
    genome = organism_id,
    trna_type = tRNA_type,
    anticodon = anticodon,
    
    # Standardize pseudogene call
    pseudo = ifelse(note == "pseudo", 1, 0)
  )

# --------------------------------
# Clean manuscript annotation data
# --------------------------------

manuscript_clean <- manuscript_trimmed_filtered %>%
  transmute(
    genome = unique.merged.organism_id.,
    trna_type = tRNA_Type,
    anticodon = Anticodon,
    
    # Standardize pseudogene call
    pseudo = ifelse(is.na(Pseudogene), 0, 1)
  )
#---------------------------------
# ----------------------------------------
# Genomes manually curated in manuscript
# ----------------------------------------
manual_genomes <- unique(manuscript_clean$genome)
#write_csv(as.data.frame(manual_genomes),"Data/manuscript_match_ids.csv")
# ----------------------------------------
# Remove original tRNAs for those genomes
# ----------------------------------------
trna_without_manual <- trna %>%
  filter(!organism_id %in% manual_genomes) %>%
  mutate(source = "trnascan")

# ----------------------------------------
# Convert manuscript rows into trna format
# ----------------------------------------
manuscript_as_trna <- manuscript_clean %>%
  transmute(
    organism = paste0(genome, "_MANUAL_CURATED"),
    tRNA_type = trna_type,
    anticodon = anticodon,
    
    # High score so they survive filtering
    inf_score = 999,
    
    note = ifelse(
      pseudo == 1,
      "pseudo",
      "manual_curated"
    ),
    
    organism_id = genome,
    
    # provenance
    source = "manuscript"
  )

# ----------------------------------------
# Combine back together
# ----------------------------------------
trna_updated <- bind_rows(
  trna_without_manual,
  manuscript_as_trna
)

#write_csv(trna_updated,"Data/trna_summary_u.csv")

#---------------------------------

# Your predictions
trna_counts <- trna_clean %>%
  count(
    genome,
    trna_type,
    anticodon,
    pseudo,
    name = "trnascan_count"
  )

# Manuscript annotations
manuscript_counts <- manuscript_clean %>%
  count(
    genome,
    trna_type,
    anticodon,
    pseudo,
    name = "manuscript_count"
  )

comparison <- full_join(
  trna_counts,
  manuscript_counts,
  by = c(
    "genome",
    "trna_type",
    "anticodon",
    "pseudo"
  )
) %>%
  mutate(
    trnascan_count = replace_na(trnascan_count, 0),
    manuscript_count = replace_na(manuscript_count, 0),
    
    difference = trnascan_count - manuscript_count
  )

perfect_matches <- comparison %>%
  filter(difference == 0)

missing_in_trnascan <- comparison %>%
  filter(difference < 0)

additional_in_trnascan <- comparison %>%
  filter(difference > 0)

genome_summary <- comparison %>%
  group_by(genome) %>%
  summarise(
    total_difference = sum(abs(difference)),
    perfect = all(difference == 0),
    n_disagreements = sum(difference != 0)
  ) %>%
  arrange(desc(n_disagreements))




missing_tbl <- missing_in_trnascan %>%
  filter(manuscript_count > 0)

extra_tbl <- additional_in_trnascan %>%
  filter(trnascan_count > 0)

candidate_pairs <- missing_tbl %>%
  inner_join(
    extra_tbl,
    by = "genome",
    suffix = c("_missing", "_extra")
  )

candidate_pairs <- candidate_pairs %>%
  mutate(
    
    same_type = trna_type_missing == trna_type_extra,
    same_anticodon = anticodon_missing == anticodon_extra,
    same_pseudo = pseudo_missing == pseudo_extra,
    
    reconciliation_type = case_when(
      
      # pseudogene disagreement
      same_type & same_anticodon & !same_pseudo ~
        "pseudo_disagreement",
      
      # anticodon mismatch
      same_type & !same_anticodon ~
        "anticodon_mismatch",
      
      # isotype mismatch
      !same_type & same_anticodon ~
        "isotype_mismatch",
      
      TRUE ~
        "other"
    )
  )

likely_same_trna <- candidate_pairs %>%
  filter(reconciliation_type != "other")

write_csv(likely_same_trna,"Data/manuscript_validation/potential_misreads_trnascan.csv")
write_csv(missing_in_trnascan,"Data/manuscript_validation/missing_in_trnascan.csv")
write_csv(additional_in_trnascan,"Data/manuscript_validation/additional_in_trnascan.csv")
