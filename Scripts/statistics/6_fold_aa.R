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

#-----------------------------------------------------------------
# Check which codons have sufficient variation and filter outliers
#-----------------------------------------------------------------

table(six_fold_subsets$Leu$codon, six_fold_subsets$Leu$CM) # CUG & UUG vary with 45 and 40 samples
six_fold_subsets$Leu <- six_fold_subsets$Leu %>%           # filter super wobble cases and CUU matches
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "CUU" & CM == "M" )
  ) %>%
  ungroup()
table(six_fold_subsets$Ser$codon, six_fold_subsets$Ser$CM) # UCG varies with 36 samples
six_fold_subsets$Ser <- six_fold_subsets$Ser %>%           # filter super wobble cases and GUU matches
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "UCU" & CM == "M") &
      !any(codon == "AGU" & CM == "M")
  ) %>%
  ungroup()
table(six_fold_subsets$Arg$codon, six_fold_subsets$Arg$CM) # AGG & CGG vary with 22 and 167 samples
six_fold_subsets$Arg <- six_fold_subsets$Arg %>%           # filter super wobble cases and GUU matches
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "CGU" & CM == "M")
  ) %>%
  ungroup()

#-----------------------------------------------------------------
# split organisms based on the state of their varying codons
#-----------------------------------------------------------------
six_fold_states <- list()

six_fold_states[["Ser"]] <- fourfold_prep(six_fold_subsets$Ser,"UCG")

leu_state1 <- six_fold_subsets$Leu %>%
  filter(codon == "CUG") %>%
  select(organism_id,
         CUG_state = CM)

leu_state2 <- six_fold_subsets$Leu %>%
  filter(codon == "UUG") %>%
  select(organism_id,
         UUG_state = CM)


leu_states <- leu_state1 %>%
  left_join(
    leu_state2,
    by = "organism_id"
  ) 

leu_states %>%
  count(CUG_state, UUG_state)

leu_states <- leu_states %>%
  filter(
  (CUG_state == "GU" & UUG_state == "GU") |
    (CUG_state == "M"  & UUG_state == "M")
)

arg_state1 <- six_fold_subsets$Arg %>%
  filter(codon == "AGG") %>%
  select(organism_id,
         AGG_state = CM)

arg_state2 <- six_fold_subsets$Arg %>%
  filter(codon == "CGG") %>%
  select(organism_id,
         CGG_state = CM)


arg_states <- arg_state1 %>%
  left_join(
    arg_state2,
    by = "organism_id"
  ) 

arg_states %>%
  count(AGG_state, CGG_state)

six_fold_states[["Leu"]] <- six_fold_subsets$Leu %>%
  inner_join(leu_states, by = "organism_id") %>%
  select(
    organism_id,
    codon,
    RSCU,
    GC3,
    CUG_state
  ) 

colnames(six_fold_states$Leu)[5] <- "state"

six_fold_states[["Arg"]] <- six_fold_subsets$Arg %>%
  inner_join(arg_states, by = "organism_id") %>%
  select(
    organism_id,
    codon,
    RSCU,
    GC3,
    AGG_state,
    CGG_state
  ) 

six_fold_states$Arg <- six_fold_states$Arg %>%
  mutate(
    state = case_when(
      CGG_state == "M" ~ "M",
      CGG_state == "GU" & AGG_state == "M" ~ "GU",
      CGG_state == "GU" & AGG_state == "GU" ~ "GU_GU"
    ),
    state = factor(state,levels = c("GU_GU","GU","M"))
  ) 

six_fold_states$Arg
six_fold_states$Leu
six_fold_states$Ser


state_rscu_boxplot(six_fold_states$Ser, "Ser")
state_rscu_boxplot(six_fold_states$Leu, "Leu")
state_rscu_boxplot(six_fold_states$Arg, "Arg")

for (i in 1:length(six_fold_states)) {
  six_fold_states[[i]] <- six_fold_states[[i]] %>%
    pivot_wider(
      names_from = codon,
      values_from = RSCU
    )
}

#-----------------------------------------------------------------
# split respons variable in 2 axes, CG vs AU; and variable codon vs the rest
#-----------------------------------------------------------------
six_fold_ready <- list()

six_fold_ready[["Ser"]] <- six_fold_axes(six_fold_states$Ser,"UCG")
six_fold_ready[["Leu"]] <- six_fold_axes(six_fold_states$Leu,c("CUG","UUG"))
six_fold_ready[["Arg"]] <- six_fold_axes(six_fold_states$Arg,c("AGG","CGG"))


for (i in 1:length(six_fold_ready)) {
  rownames(six_fold_ready[[i]]) <- six_fold_ready[[i]]$organism_id
}

#-----------------------------------------------------------------
# fit the models
#-----------------------------------------------------------------

# ser
fit_ser <- analyze_sixfold(six_fold_ready$Ser,tree)
best_ser <- best_sixfold(fit_ser)
best_ser$gc_delta <-
  delta_aic(best_ser$gc_aic)

best_ser$var_delta <-
  delta_aic(best_ser$var_aic)

# Leu
fit_leu <- analyze_sixfold(six_fold_ready$Leu,tree)
best_leu <- best_sixfold(fit_leu)
best_leu$gc_delta <-
  delta_aic(best_leu$gc_aic)

best_leu$var_delta <-
  delta_aic(best_leu$var_aic)

# Arg
fit_arg <- analyze_sixfold(six_fold_ready$Arg,tree)
best_arg <- best_sixfold(fit_arg)
best_arg$gc_delta <-
  delta_aic(best_arg$gc_aic)

best_arg$var_delta <-
  delta_aic(best_arg$var_aic)

# save output
save_sixfold(
  fit = fit_ser,
  best = best_ser,
  aa = "Plots/aa/six_fold/Ser"
)

save_sixfold(
  fit = fit_leu,
  best = best_leu,
  aa = "Plots/aa/six_fold/Leu"
)

save_sixfold(
  fit = fit_arg,
  best = best_arg,
  aa = "Plots/aa/six_fold/Arg"
)

