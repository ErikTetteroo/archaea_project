#load packages
source("R/load_project.R")

#load in codon pair dataset & merged dataset
codon_pair_usage <- read_csv("Data/Raw_data/combined_codon_pair_usage.csv")
merged <- read_csv("Data/cleaned_data/merged_codon_usage_data.csv")

#standardize codon format
merged <- merged %>%
  mutate(codon = gsub("U", "T", codon))

#create lookup table
cv_lookup <- merged %>%
  select(organism_id, codon, CV)

#add coverage value codon 1
codon_pair_usage <- codon_pair_usage %>%
  left_join(
    cv_lookup,
    by = c("organism" = "organism_id", "c1" = "codon")
  ) %>%
  rename(CV1 = CV)

#add coverage value codon 2
codon_pair_usage <- codon_pair_usage %>%
  left_join(
    cv_lookup,
    by = c("organism" = "organism_id", "c2" = "codon")
  ) %>%
  rename(CV2 = CV)

#combine coverage values (currently multiplied)
codon_pair_usage <- codon_pair_usage %>%
  mutate(CVP = CV1 * CV2)

#take log of both
codon_pair_usage <- codon_pair_usage %>%
  mutate(log2_CVP = log2(CVP),
         log2_RDCU = log2(RDCU))

#plot RDCU vs coverage value (Takes long and is not very insightful atm)
#ggplot(Codon_pair_usage, aes(x = log2_CVP, y = log2_RDCU)) +
#  geom_point(alpha = 0.1, size = 0.5) +
#  theme_minimal() +
#  labs(
#    title = "Relationship between codon pair coverage and RDCU",
#    x = "log2(CVP)",
#    y = "log2(RDCU)"
#  )

#many 0's in coverage value stop further analysis atm



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
