n_total <- length(unique(merged$organism_id))
n_in_tree <- nrow(present_in_tree)
n_missing <- nrow(missing_in_tree)

coverage <- n_in_tree / n_total
coverage

missing_by_phylum <- missing_in_tree %>%
  count(phylum, sort = TRUE)

present_by_phylum <- present_in_tree %>%
  count(phylum, sort = TRUE)

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

merged %>%
  mutate(in_tree = organism_id %in% present_in_tree$organism_id) %>%
  group_by(in_tree) %>%
  summarise(
    mean_GC3 = mean(GC3, na.rm = TRUE),
    sd_GC3 = sd(GC3, na.rm = TRUE),
    mean_total = mean(total, na.rm = TRUE)
  )

merged2 <- rbind(present_in_tree,missing_in_tree)

tree_phyla <- data.frame(
  organism_id = tree$tip.label
) %>%
  left_join(merged2, by = "organism_id")

table(tree_phyla$phylum)

table(merged2$phylum)

phylum_loss <- merged2 %>%
  count(phylum) %>%
  rename(total = n) %>%
  left_join(tree_phyla %>% count(phylum) %>% rename(in_tree = n),
            by = "phylum") %>%
  mutate(in_tree = ifelse(is.na(in_tree), 0, in_tree),
         loss_fraction = 1 - in_tree / total)
