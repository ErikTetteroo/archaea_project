# load project
source("R/load_project.R")

# read in files
merged <- read_csv("Data/cleaned_data/codon_usage_data_stats_ready.csv")
tree <- read.tree("Data/cleaned_data/cleaned_archaea.tree")
valid_amino_acids <- read_csv("Data/cleaned_data/aa_to_analyse.csv")

# subset the interesting 4 fold amino acids
six_fold <- valid_amino_acids$aa[valid_amino_acids$Freq==6]

six_fold_subsets <- list()

for (n in six_fold) {
  six_fold_subsets[[n]] <- aa_subset(merged, n)
}

table(six_fold_subsets$Leu$codon, six_fold_subsets$Leu$CM)
six_fold_subsets$Leu <- six_fold_subsets$Leu %>%           # filter super wobble cases and CUU matches
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "CUU" & CM == "M" )
  ) %>%
  ungroup()
table(six_fold_subsets$Ser$codon, six_fold_subsets$Ser$CM)
six_fold_subsets$Ser <- six_fold_subsets$Ser %>%           # filter super wobble cases and GUU matches
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "UCU" & CM == "M") &
      !any(codon == "AGU" & CM == "M")
  ) %>%
  ungroup()
table(six_fold_subsets$Arg$codon, six_fold_subsets$Arg$CM)
six_fold_subsets$Arg <- six_fold_subsets$Arg %>%           # filter super wobble cases and GUU matches
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "CGU" & CM == "M")
  ) %>%
  ungroup()
