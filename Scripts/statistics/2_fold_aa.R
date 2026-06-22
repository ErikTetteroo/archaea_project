# load project
source("R/load_project.R")

# read in files
merged <- read_csv("Data/cleaned_data/codon_usage_data_stats_ready.csv")
tree <- read.tree("Data/cleaned_data/cleaned_archaea.tree")
valid_amino_acids <- read_csv("Data/cleaned_data/aa_to_analyse.csv")

# subset the interesting 2 fold amino acids
two_fold <- valid_amino_acids$aa[valid_amino_acids$Freq==2]

two_fold_subsets <- list()

for (n in two_fold) {
  two_fold_subsets[[n]] <- aa_subset(merged, n)
}


# check the different codon coverage's
table(two_fold_subsets$Gln$codon, two_fold_subsets$Gln$CM)  # CAG varies with 74 alternative samples
table(two_fold_subsets$Lys$codon, two_fold_subsets$Lys$CM)  # AAG varies with 41 alternative samples
table(two_fold_subsets$Glu$codon, two_fold_subsets$Glu$CM)  # GAG varies with 42 alternative samples

# fit models
gln <- analyze_twofold(two_fold_subsets$Gln,"CAG",tree)
lys <- analyze_twofold(two_fold_subsets$Lys,"AAG",tree)
glu <- analyze_twofold(two_fold_subsets$Glu,"GAG",tree)

# save output
save_twofold(gln,"plots/aa/Gln")
save_twofold(lys,"plots/aa/Lys")
save_twofold(glu,"plots/aa/Glu")

# summary table
bind_rows(
  summary_table(gln, "Gln"),
  summary_table(lys, "Lys"),
  summary_table(glu, "Glu")
)
