#-------------------------------------------------------------------------------
# manuscript validation
#-------------------------------------------------------------------------------

# step 1: run this script to check which ids match with our data
source('Scripts/data_prep/manuscript/match_to_manuscript.R')

# step 2: run the manuscript_trna_update.R script to overwrite trnascan results 
#         with the more accurate trna data from the manuscript 
#         (for the matching ids)
source('Scripts/data_prep/manuscript/manuscript_trna_update.R')

#-------------------------------------------------------------------------------
# merge trna and codon usage datasets
#-------------------------------------------------------------------------------

# step 1: run this script to create trna coverage table and plot
source('Scripts/data_prep/tRNA_coverage.R')

# step 2: run the data_prep/Merging.R script to create merged 
#         codon_usage/coverage dataframe
source("Scripts/data_prep/Merging.R")

# deal with illegal 'no coverage' codons

#-------------------------------------------------------------------------------
# deviations from standard set approach 
# (should result in same outcome as steps above)
#-------------------------------------------------------------------------------

#tRNA_dive.R
# step 1: check/count rare trnas and missing common trnas
source("Scripts/data_prep/tRNA_dive.R")

# step 2: see what the coverage is of the 46 tRNA standard set of archaea
source("Scripts/data_prep/tRNA_coverage_of_standard_set.R")

# step 3: map out the deviations from the standard set for organisms
source("Scripts/data_prep/deviations_of_standard_set.R")

# step 4: plot the codon coverage of these deviated trna pools to illustrate 
#         impossibilities (codons with no coverage at all)
source("Scripts/creating_plots/coverage_plot_deviations.R")

# step 5: compare coverage of this deviation approach with earlier merge attempt
source("Scripts/data_prep/deviation_vs_direct_validation.R")

#-------------------------------------------------------------------------------
# taxonomy
#-------------------------------------------------------------------------------

# This script was used to clean tree tip labels and to explore which organisms
# were cut in the tree. Originally a pcoa was used to add tree data to the 
# merged dataset. But this didn't provide enought taxonomic coverage for the 
# analysis and the code for this was therefore moved to the outdated scripts
source("Scripts/data_prep/taxonomy.R")

#-------------------------------------------------------------------------------
# statistics
#-------------------------------------------------------------------------------

# step 1: Check which amino acids are viable for analysis and visualize this
source("Scripts/statistics/stats_prep.R")

# step 2: analyse the twofold amino acids
source("Scripts/statistics/2_fold_aa.R")

# step 3: analyse the fourfold amino acids
source("Scripts/statistics/4_fold_aa.R")

# step 4: analyse the sixfold amino acids
source("Scripts/statistics/6_fold_aa.R")

#-------------------------------------------------------------------------------
# codon pairs
#-------------------------------------------------------------------------------

