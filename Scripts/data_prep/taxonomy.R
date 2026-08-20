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


plot(
  tree,
  cex = 0.4
)

tiff(
  "Plots/phylogenetic_tree.tiff",
  width = 2000,
  height = 2800,
  res = 300
)

plot(
  tree,
  show.tip.label = FALSE
)

dev.off()

tax <- merged[merged$in_tree,]

library(dplyr)

tax_org <- tax %>%
  select(
    organism_id,
    organism,
    phylum,
    class,
    order,
    family,
    genus
  ) %>%
  distinct()

tax_org %>%
  summarise(
    organisms = n_distinct(organism_id),
    phyla = n_distinct(phylum),
    classes = n_distinct(class),
    orders = n_distinct(order),
    families = n_distinct(family),
    genera = n_distinct(genus)
  )

tax_summary <- bind_rows(
  
  tax_org %>%
    count(phylum, name = "n") %>%
    mutate(
      level = "Phylum",
      taxon = phylum
    ) %>%
    select(level, taxon, n),
  
  tax_org %>%
    count(class, name = "n") %>%
    mutate(
      level = "Class",
      taxon = class
    ) %>%
    select(level, taxon, n),
  
  tax_org %>%
    count(order, name = "n") %>%
    mutate(
      level = "Order",
      taxon = order
    ) %>%
    select(level, taxon, n),
  
  tax_org %>%
    count(family, name = "n") %>%
    mutate(
      level = "Family",
      taxon = family
    ) %>%
    select(level, taxon, n),
  
  tax_org %>%
    count(genus, name = "n") %>%
    mutate(
      level = "Genus",
      taxon = genus
    ) %>%
    select(level, taxon, n)
  
) %>%
  group_by(level) %>%
  mutate(
    percentage = 100 * n / sum(n)
  ) %>%
  ungroup() %>%
  arrange(
    factor(level, levels = c(
      "Phylum", "Class", "Order", "Family", "Genus"
    )),
    desc(n)
  )

tax_summary
