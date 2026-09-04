merged <- read_csv("Data/Merged_data_m.csv")
mcc <- read_csv("Data/missing_common_combos.csv")
prc <- read_csv("Data/present_rare_combos.csv")

mergedo <- merged[order(merged$GC3),]
unique(merged$species)

head(prc)

library(dplyr)

# organisms present in both datasets
shared_orgs <- intersect(
  unique(mcc$organism_id),
  unique(prc$organism_id)
)

# missing combos
mcc_subset <- mcc %>%
  filter(organism_id %in% shared_orgs) %>%
  transmute(
    tRNA_type,
    anticodon,
    present = FALSE,
    organism_id,
    notes
  )

# present rare combos
prc_subset <- prc %>%
  filter(organism_id %in% shared_orgs) %>%
  transmute(
    tRNA_type,
    anticodon,
    present = TRUE,
    organism_id,
    notes = NA_character_
  )

combined_trna_changes <- bind_rows(
  mcc_subset,
  prc_subset
) %>%
  arrange(organism_id, present, tRNA_type, anticodon)

unique(prc_subset$organism_id)

combined_trna_changes[combined_trna_changes$organism_id==unique(combined_trna_changes$organism_id)[44],]
