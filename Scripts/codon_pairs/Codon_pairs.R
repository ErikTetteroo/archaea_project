#load packages
source("R/load_project.R")

#load in codon pair dataset & merged dataset
codon_pair_usage <- read_csv("Data/Raw_data/combined_codon_pair_usage.csv")
merged <- read_csv("Data/cleaned_data/merged_codon_usage_data.csv")

#standardize codon format
merged$codon <- rna_to_dna(merged$codon)

#create lookup table
cm_lookup <- merged %>%
  select(organism_id, codon, CM)


gc_lookup <- merged %>%
  select(organism_id, codon, GC3)


#add GC3
codon_pair_usage <- codon_pair_usage %>%
  left_join(
    gc_lookup,
    by = c("organism" = "organism_id", "c1" = "codon")
  )

#add coverage value codon 1
codon_pair_usage <- codon_pair_usage %>%
  left_join(
    cm_lookup,
    by = c("organism" = "organism_id", "c1" = "codon")
  ) %>%
  rename(CM1 = CM)

#add coverage value codon 2
codon_pair_usage <- codon_pair_usage %>%
  left_join(
    cm_lookup,
    by = c("organism" = "organism_id", "c2" = "codon")
  ) %>%
  rename(CM2 = CM)

#combine coverage values (currently multiplied)
codon_pair_usage <- codon_pair_usage %>%
  mutate(CMP = paste0(CM1,"_X_",CM2))

#take the log of RDCU
codon_pair_usage <- codon_pair_usage %>%
  mutate(log2_RDCU = log2(RDCU))

order <- sort(table(codon_pair_usage$CMP), decreasing = TRUE)

labs <- paste0(names(order), " (n=", order, ")")
names(labs) <- names(order)

ggplot(
  codon_pair_usage,
  aes(
    x = factor(CMP, levels = names(order)),
    y = log2_RDCU
  )
) +
  geom_boxplot() +
  scale_x_discrete(labels = labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

write_csv(codon_pair_usage,"Data/raw_data/codon_pair_usage_clean.csv" )

#creating a heatmap#######################################################

#extract relevant columns
Codon_pair_usage_hm <- codon_pair_usage[,c(1,2,5)]

#filter interesting pairs by log distance of RDCU 1
top_pairs <- Codon_pair_usage_hm %>%
  group_by(codon_pair) %>%
  summarise(score = mean(abs(log2(RDCU)), na.rm = TRUE)) %>%
  arrange(desc(score)) %>%
  slice_head(n = 200) %>%     #decide how many pairs you'd like to include
  pull(codon_pair)

df_small <- Codon_pair_usage_hm %>%
  filter(codon_pair %in% top_pairs)

#create matrix
matrix_df <- df_small %>%
  pivot_wider(
    names_from = codon_pair,
    values_from = RDCU
  )

#log transform
matrix_log <- matrix_df %>%
  mutate(across(-organism, ~log2(ifelse(. <= 0, NA, .))))

#prep for plotting
plot_df <- matrix_log %>%
  pivot_longer(
    -organism,
    names_to = "codon_pair",
    values_to = "RDCU"
  )


#order rows/columns by similarity
mat <- matrix_log %>%
  column_to_rownames("organism") %>%
  as.matrix()

row_order <- hclust(dist(mat))$order
ordered_organisms <- rownames(mat)[row_order]

col_order <- hclust(dist(t(mat)))$order
ordered_pairs <- colnames(mat)[col_order]

plot_df$organism <- factor(plot_df$organism, levels = ordered_organisms)
plot_df$codon_pair <- factor(plot_df$codon_pair, levels = ordered_pairs)

#plot heatmap
ggplot(plot_df, aes(x = codon_pair, y = organism, fill = RDCU)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # hide (too many)
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 5)
  ) +
  labs(
    title = "RDCU Heatmap (Top 500 Codon Pairs, log2 scale)",
    x = "Codon Pair",
    y = "Organism",
    fill = "log2(RDCU)"
  )
