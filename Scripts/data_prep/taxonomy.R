#load project
source("R/load_project.R")

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

missing_in_tree <- distinct(merged[merged$organism_id %in% 
                                     unique(merged$organism_id[!merged$in_tree]),c(1,17:23)])
present_in_tree <- distinct(merged[merged$organism_id %in% 
                                     unique(merged$organism_id[merged$in_tree]),c(1,17:23)])

# write cleaned tree
write.tree(tree, "Data/cleaned_data/cleaned_archaea.tree")
write_csv(merged[merged$in_tree,], "Data/cleaned_data/codon_usage_data_stats_ready.csv")


#---------------------------------------------------------------------
# taxonomy summary stats
#---------------------------------------------------------------------

n_total <- length(unique(merged$organism_id))
n_in_tree <- nrow(present_in_tree)
n_missing <- nrow(missing_in_tree)

# percentage of organisms covered
coverage <- n_in_tree / n_total
coverage

missing_by_phylum <- missing_in_tree %>%
  count(phylum, sort = TRUE)

present_by_phylum <- present_in_tree %>%
  count(phylum, sort = TRUE)

# coverage per phylum
coverage_phylum <- full_join(
  present_by_phylum %>% rename(present = n),
  missing_by_phylum %>% rename(missing = n),
  by = "phylum"
) %>%
  mutate(
    present = ifelse(is.na(present), 0, present),
    missing = ifelse(is.na(missing), 0, missing),
    coverage = present / (present + missing)
  ) %>%
  arrange(coverage)
coverage_phylum

# metrics of organisms cut from tree
merged %>%
  mutate(in_tree = organism_id %in% present_in_tree$organism_id) %>%
  group_by(in_tree) %>%
  summarise(
    mean_GC3 = mean(GC3, na.rm = TRUE),
    sd_GC3 = sd(GC3, na.rm = TRUE),
    mean_total = mean(total, na.rm = TRUE)
  )

