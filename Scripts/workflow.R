#-------
# manuscript validation
#-------

# step 1: run the match_to_manuscript.R script to check which ids match between datasets
source('Scripts/data_prep/manuscript/match_to_manuscript.R')

# step 2: run the manuscript_trna_update.R script to overwrite trnascan results with the 
#         more accurate trna data from the manuscript (for the matching ids)
source('Scripts/data_prep/manuscript/manuscript_trna_update.R')

#-------
# merge trna and codon usage datasets
#-------

# step 1: run this script to create trna coverage table and plot
source('Scripts/data_prep/tRNA_coverage.R')

# step 2: run the data_prep/Merging.R script to create merged codon_usage/coverage dataframe
source("Scripts/data_prep/Merging.R")

# deal with illegal 'no coverage' codons

#-------
# deviations from standard set approach (should result in same outcome as steps above)
#-------

#tRNA_dive.R
# Step 1: check/count rare trnas and missing common trnas
source("Scripts/data_prep/tRNA_dive.R")

# Step 2: see what the coverage is of the 46 tRNA standard set of archaea
source("Scripts/data_prep/tRNA_coverage_of_standard_set.R")

#deviations_of_standard_set.R
source("Scripts/data_prep/deviations_of_standard_set.R")

#-------
# codon pairs
#-------

