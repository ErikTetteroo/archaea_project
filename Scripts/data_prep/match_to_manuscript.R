source("R/load_project.R")

#load in data
merged <- read_csv("Data\\Merged_data.csv")
acc <- read_csv("Data\\Raw_data\\Assembly_Acc_manuscript.txt")

#check which genomes were also used in manuscript
u_acc <- unique(acc)
u_acc <- u_acc %>%
  mutate(core = sub("^GCA_", "", Assembly_Acc))

merged_unique <- as.data.frame(unique(merged$organism_id)) %>%
  mutate(core = sub("^GCF_", "", unique(merged$organism_id)))

matches <- inner_join(u_acc, merged_unique, by = "core")

#save matches lookup table and add in a column to the merged dataframe
write_csv(matches,file = "Data\\manuscript_match_id.csv")
merged$in_mnscrpt <- ifelse(merged$organism_id %in% matches$`unique(merged$organism_id)`, 1, 0)
write_csv(merged,file = "Data\\Merged_data_m.csv")
