#load project
source("R/load_project.R")

# read in data
merged <- read_csv("Data/cleaned_data/codon_usage_data_stats_ready.csv")
tree <- read.tree("Data/cleaned_data/cleaned_archaea.tree")

# add a phylogeney proxy through pcoa
dist <- cophenetic(tree)

pcoa <-pcoa(dist)

phy <- data.frame(
  organism_id = rownames(pcoa$vectors),
  phyPC1 = pcoa$vectors[,1],
  phyPC2 = pcoa$vectors[,2],
  phyPC3 = pcoa$vectors[,3],
  phyPC4 = pcoa$vectors[,4]
)

merged <- merged %>%
  left_join(phy, by = "organism_id")

# output updated merged file
write_csv(merged,"Data/cleaned_data/merged_codon_usage_data_with_pcoa_taxonomy.csv")
