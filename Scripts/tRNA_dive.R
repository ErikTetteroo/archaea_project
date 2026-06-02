library(dplyr)
library(tidyverse)
library(tidyr)

#load in data
tRNA <- read_csv("Data\\trna_summary_u.csv")

# Final filtered dataset
trna_filtered <- tRNA %>%
  filter(
    inf_score > 30,
    is.na(note) | !str_detect(note, "pseudo"),
    !str_detect(anticodon, "N")
  )

# Everything filtered out
trna_low_quality <- tRNA %>%
  filter(
    inf_score <= 30 |
      (!is.na(note) & str_detect(note, "pseudo")) |
      str_detect(anticodon, "N")
  )

trna_clean <- trna_filtered[,c(6,2,3)]

#How many trnas does each organism contain total
total_trna_per_org <- trna_clean %>%
  group_by(organism_id) %>%
  summarise(total_trna = n(), .groups = "drop")

#How many unique trnas does each organism contain
unique_trna_per_org <- trna_clean %>%
  distinct(organism_id, tRNA_type, anticodon) %>%
  group_by(organism_id) %>%
  summarise(unique_trna = n(), .groups = "drop")

#Combine metrics
summary_table <- total_trna_per_org %>%
  left_join(unique_trna_per_org, by = "organism_id")


#count how often each trna type is present in each organism
trna_combo_counts <- trna_clean %>%
  distinct(organism_id, tRNA_type, anticodon) %>%  # remove within-organism duplicates
  group_by(tRNA_type, anticodon) %>%
  summarise(n_organisms = n(), .groups = "drop") %>%
  arrange(desc(n_organisms))


#find out which organisms have rare combos or are missing common combos
trna_presence <- trna_clean %>%
  distinct(organism_id, tRNA_type, anticodon)





rare_combos <- trna_combo_counts %>%
  filter(n_organisms <= 21)

rare_combo_organisms <- rare_combos %>%
  inner_join(trna_presence, by = c("tRNA_type", "anticodon")) %>%
  arrange(tRNA_type, anticodon)

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


#check if missing common combos can be explained by low quality trna's
missing_with_qc <- missing_common_combos %>%
  left_join(
    trna_low_quality %>%
      mutate(is_pseudo = note == "pseudo"),
    by = c("organism_id", "tRNA_type", "anticodon")
  ) %>%
  group_by(organism_id, tRNA_type, anticodon) %>%
  summarise(
    M = sum(!is.na(note)),              
    P = sum(is_pseudo, na.rm = TRUE),   
    .groups = "drop"
  )

n_detection <- trna_low_quality %>%
  mutate(has_N = grepl("N", anticodon)) %>%
  group_by(organism_id) %>%
  summarise(
    N_hits = sum(has_N),  
    N_unique = n_distinct(anticodon[has_N]),
    .groups = "drop"
  ) %>%
  filter(N_hits > 0)

missing_with_qc <- missing_with_qc %>%
  left_join(
    n_detection %>% select(organism_id, N_hits),
    by = "organism_id"
  ) %>%
  mutate(
    N_hits = coalesce(N_hits, 0)
  )

missing_with_qc$U <- missing_with_qc$M + missing_with_qc$N_hits

#manually check cases
missing_check <- missing_with_qc[missing_with_qc$U>0,]
missing_check$notes <- "-"

###Long and painful process!!
#source("Functions\\manual_check.R")
#orgs <- unique(missing_check$organism_id)
#for(i in seq_along(orgs)) {
#
#  review_missing(orgs[i])
#
#  readline(prompt = "Press [enter] for next organism...")
#}


#attach the manually written notes back onto the data
lines <- readLines("missing_common_combos_manual_check.txt")

lines <- trimws(lines)
lines <- lines[lines != ""]

library(stringr)

#parse notes
parsed_entries <- list()

i <- 1

while(i <= length(lines)) {
  
  # organism_id line
  org_id <- lines[i]
  i <- i + 1
  
  # safety check
  if(i > length(lines)) break
  
  # next line determines entry type
  next_line <- lines[i]
  
  # CASE 1:
  # organism + note only
  if(str_detect(next_line, '^".*"$')) {
    
    note <- str_remove_all(next_line, '^"|"$')
    
    parsed_entries[[length(parsed_entries) + 1]] <- tibble(
      organism_id = org_id,
      tRNA_type = NA_character_,
      anticodon = NA_character_,
      note = note
    )
    
    i <- i + 1
    
  } else {
    
    # CASE 2:
    # organism + tRNA + anticodon + note
    
    trna_type <- next_line
    i <- i + 1
    
    if(i > length(lines)) break
    anticodon <- lines[i]
    i <- i + 1
    
    if(i > length(lines)) break
    note <- str_remove_all(lines[i], '^"|"$')
    i <- i + 1
    
    parsed_entries[[length(parsed_entries) + 1]] <- tibble(
      organism_id = org_id,
      tRNA_type = trna_type,
      anticodon = anticodon,
      note = note
    )
  }
}

manual_notes <- bind_rows(parsed_entries)

#attach organism notes
organism_notes <- manual_notes %>%
  filter(is.na(tRNA_type))

missing_check2 <- missing_check %>%
  left_join(
    organism_notes %>%
      select(organism_id, organism_note = note),
    by = "organism_id"
  ) %>%
  mutate(
    notes = ifelse(
      is.na(organism_note),
      notes,
      organism_note
    )
  ) %>%
  select(-organism_note)

#attach notes for specific combos
combo_notes <- manual_notes %>%
  filter(!is.na(tRNA_type))

missing_check2 <- missing_check2 %>%
  left_join(
    combo_notes %>%
      rename(combo_note = note),
    by = c("organism_id", "tRNA_type", "anticodon")
  ) %>%
  mutate(
    notes = ifelse(
      is.na(combo_note),
      notes,
      combo_note
    )
  ) %>%
  select(-combo_note)

#make a copy of parsed manually written notes
#write.csv(manual_notes,
#          "parsed_manual_notes.csv",
#          row.names = FALSE)

#attach notes to larger dataframe
missing_common_combos_annotated <- missing_common_combos %>%
  left_join(
    missing_check2 %>%
      select(organism_id, tRNA_type, anticodon, notes),
    by = c("organism_id", "tRNA_type", "anticodon")
  )

#save a copy of both
write_csv(missing_common_combos,"Data\\missing_common_combos_u.csv")
write_csv(rare_combo_organisms,"Data\\present_rare_combos_u.csv")
