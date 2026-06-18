#load project
source("R/load_project.R")
library(ape)

# read in data
merged <- read_csv("Data/cleaned_data/merged_codon_usage_data.csv")
tree <- read.tree("Data/tree/filtered_690.tree")


# reformat tree tip labels
tree$tip.label <- sub("^RS_", "", tree$tip.label)


# check how many organisms are missing from the list
length(tree$tip.label)
length(unique(merged$organism_id))

sum(unique(merged$organism_id) %in% tree$tip.label)


#prune merged down to only ids present in tree
common_ids <- intersect(merged$organism_id, tree$tip.label)
merged_tree <- merged[merged$organism_id %in% common_ids, ]

merged$in_tree <- merged$organism_id %in% tree$tip.label

missing_in_tree <- distinct(merged[merged$organism_id %in% unique(merged$organism_id[!merged$in_tree]),c(1,17:23)])
present_in_tree <- distinct(merged[merged$organism_id %in% unique(merged$organism_id[merged$in_tree]),c(1,17:23)])

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
