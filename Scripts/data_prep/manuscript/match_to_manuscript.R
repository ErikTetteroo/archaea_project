local({
  source("R/load_project.R")

#load in data
my_data <- read_csv("Data/Raw_data/combined_codon_usage.csv")
manuscript_accesions <- read_csv("Data/Raw_data/Assembly_Acc_manuscript.txt")

#check which genomes were also used in manuscript
u_accesions <- unique(manuscript_accesions)
u_accesions <- u_accesions %>%
  mutate(core = sub("^GCA_", "", Assembly_Acc))

organism_id <- unique(my_data$organism)
my_data_unique <- as.data.frame(organism_id) %>%
  mutate(core = sub("^GCF_", "", organism_id))

matches <- inner_join(u_accesions, my_data_unique, by = "core")

# Save lookup table of genomes shared with the manuscript
write_csv(matches,file = "Data/manuscript/manuscript_id_match.csv")
})