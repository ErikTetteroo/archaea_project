#load project
source("R/load_project.R")


# Load data
merged <- read_csv("Data/cleaned_data/codon_usage_data_stats_ready.csv")
tree <- read.tree("Data/cleaned_data/cleaned_archaea.tree")

# Match phylum to tree tips
tax_org <- merged %>%
  filter(in_tree) %>%
  select(organism_id, phylum) %>%
  distinct()

tip_phylum <- tax_org$phylum[
  match(tree$tip.label, tax_org$organism_id)
]

# Define phylum colours
phylum_cols <- c(
  Methanobacteriota = "#1B4F72",
  Thermoproteota = "#D35400",
  Nitrososphaerota = "#27AE60",
  Thermoplasmatota = "#8E44AD",
  `Candidatus Nanohalarchaeota` = "#16A085",
  Nanobdellota = "#C0392B",
  Microcaldota = "#7F8C8D",
  Promethearchaeota = "#C49A00"
)

tip_cols <- unname(phylum_cols[tip_phylum])

# Create phylogenetic tree plot

png(
  filename = "Plots/phylogenetic_tree_phylum.png",
  width = 3000,
  height = 4000,
  res = 300
)

plot(
  tree,
  show.tip.label = FALSE
)

tiplabels(
  pch = 19,
  col = tip_cols,
  cex = 0.8
)

legend(
  "topleft",
  legend = names(phylum_cols),
  col = phylum_cols,
  pch = 19,
  cex = 0.8,
  bty = "n"
)

dev.off()

# Class composition summary

class_summary <- merged %>%
  filter(in_tree) %>%
  select(organism_id, class) %>%
  distinct() %>%
  count(class, name = "n") %>%
  mutate(
    percentage = 100 * n / sum(n)
  ) %>%
  arrange(desc(n))

class_summary
