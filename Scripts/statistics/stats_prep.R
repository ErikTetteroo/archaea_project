# load project
source("R/load_project.R")

# read in files
merged <- read_csv("Data/cleaned_data/codon_usage_data_stats_ready.csv")
tree <- read.tree("Data/cleaned_data/cleaned_archaea.tree")

# identify interesting codons (with alternative coverage in sufficient numbers)
merged %>%
  group_by(codon) %>%
  summarise(
    n_M  = sum(CM == "M", na.rm = TRUE),
    n_GU = sum(CM == "GU", na.rm = TRUE),
    n_SU = sum(CM == "SU", na.rm = TRUE),
    n_M2 = sum(CM == "M2", na.rm = TRUE)
  ) %>%
  print(n=63)

# codons with at least 20 organisms with alternative coverage are deemed interesting
valid_codons <- merged %>%
  filter(!is.na(CM)) %>%
  count(codon, CM) %>%
  group_by(codon) %>%
  summarise(
    n_categories_20 = sum(n >= 20)
  ) %>%
  filter(n_categories_20 >= 2) %>%
  pull(codon)

# The amino acids that contain these codons will be analysed
table <- give_codon_table()
valid_amino_acids <- tibble(
  aa = unique(table$aa[table$codon_rna %in% valid_codons]))

valid_amino_acids <- valid_amino_acids %>%
  left_join(as.data.frame(table(table$aa)),
            by = c("aa" = "Var1"))

valid_amino_acids <- valid_amino_acids[order(valid_amino_acids$Freq),]

# plot of the interesting codons
merged_c <- merged[!is.na(merged$CM),]

cm_prop <- merged_c |>
  summarise(
    n = n(),
    .by = c(codon, CM)
  ) |>
  group_by(codon) |>
  mutate(prop = n / sum(n)) |>
  ungroup()


aa_order <- c(
  "Met", "Trp", "Asn", "Asp", "Cys",
  "Gln", "Glu", "His", "Lys", "Phe",
  "Tyr", "Ile", "Ala", "Gly", "Pro",
  "Thr", "Val", "Arg", "Leu", "Ser"
)

cm_prop2 <- cm_prop %>%
  left_join(
    table %>% select(aa, codon = codon_rna),
    by = "codon"
  ) %>%
  filter(!aa == "TER") %>%
  mutate(
    aa = factor(aa, levels = aa_order)
  )

p <- ggplot(cm_prop2, aes(codon, prop, fill = CM)) +
  geom_col() +
  facet_wrap(~ aa, scales = "free_x") +
  labs(
    x = "Codon",
    y = "Proportion",
    fill = "Coverage type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

# save plot
ggsave(
  filename = "Plots/codon_variation.png",
  plot = p,
  width = 12,
  height = 10,
  dpi = 300
)

# output amino acids to analyse
write_csv(valid_amino_acids,"Data/cleaned_data/aa_to_analyse.csv")
