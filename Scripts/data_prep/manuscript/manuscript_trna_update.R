local({source("R/load_project.R")

manuscript_data <- read_csv("Data/Raw_data/archaea_complete_curated_zenodo.csv")
trna_dat <- read.csv("Data/Raw_data/trnas_summary.csv")
matches <- read.csv("Data/manuscript/manuscript_id_match.csv")

#Create consistent organism id
trna <- trna_dat %>%
  mutate(organism_id = str_extract(organism, "GCF_[0-9]+\\.[0-9]+"))

# Add matching Assembly_Acc to trna
trna_trimmed <- trna %>%
  inner_join(
    matches,
    by = c("organism_id")
  )

# Add matching organism_id to manuscript_data
manuscript_trimmed <- manuscript_data %>%
  inner_join(
    matches,
    by = "Assembly_Acc"
  )

# filter out not_ok/maybe curations from manuscript
manuscript_trimmed_filtered <- manuscript_trimmed %>%
  filter(is.na(curation) | curation == "ok")

# manually change 1 lys to Ser as note says it is a CGA anticodon corresponding to ser
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
    genome = organism_id,
    trna_type = tRNA_Type,
    anticodon = Anticodon,
    
    # Standardize pseudogene call
    pseudo = ifelse(is.na(Pseudogene), 0, 1)
  )


# ----------------------------------------
# Remove original tRNAs for those genomes
# ----------------------------------------

trna_without_manual <- trna %>%
  filter(!organism_id %in% matches$organism_id) %>%
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

# ----------------------------------------
# correct manually identified trna mismatches
# ----------------------------------------
manual_trna_fixes <- tibble::tribble(
  ~organism_id,         ~old_type, ~old_anticodon, ~new_type, ~new_anticodon,
  "GCF_000246985.2",   "Pro",     "NGG",          "Pro",     "GGG",
  "GCF_000246985.2",   "Undet",   "NTG",          "Gln",     "CTG",
  "GCF_000246985.2",   "Leu",     "NAG",          "Leu",     "CAG"
)

trna_updated <- trna_updated %>%
  left_join(
    manual_trna_fixes,
    by = c(
      "organism_id",
      "tRNA_type" = "old_type",
      "anticodon" = "old_anticodon"
    )
  ) %>%
  mutate(
    tRNA_type = coalesce(new_type, tRNA_type),
    anticodon = coalesce(new_anticodon, anticodon)
  ) %>%
  select(-new_type, -new_anticodon)

write_csv(trna_updated,"Data/cleaned_data/trna_summary_c.csv")})
