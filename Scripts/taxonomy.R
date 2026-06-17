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
