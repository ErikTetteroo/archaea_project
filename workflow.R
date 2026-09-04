#-------------------------------------------------------------------------------
# manuscript validation and tRNA curation
#-------------------------------------------------------------------------------

# The manuscript by van de Gulik et al. contains manually curated tRNA annotations 
# for a subset of genomes. These annotations are used instead of the original 
# tRNAscan-SE predictions for matching genomes.

# Step 1: identify genomes shared between the NCBI dataset and manuscript
source("Scripts/data_prep/manuscript/match_to_manuscript.R")

# Step 2: replace tRNAscan-SE annotations with manuscript annotations for
# matching genomes and apply manual corrections
source("Scripts/data_prep/manuscript/manuscript_trna_update.R")


#-------------------------------------------------------------------------------
# merge trna and codon usage datasets
#-------------------------------------------------------------------------------

# step 1: run this script to create trna coverage table and plot
source('Scripts/data_prep/tRNA_coverage.R')

# step 2: run the data_prep/Merging.R script to create merged 
#         codon_usage/coverage dataframe
source("Scripts/data_prep/Merging.R")

#-------------------------------------------------------------------------------
# deviations from standard set approach 
# (should result in same outcome as steps above)
#-------------------------------------------------------------------------------

#tRNA_dive.R
# step 1: check/count rare trnas and missing common trnas
source("Scripts/standard_set/tRNA_dive.R")

# step 2: see what the coverage is of the 46 tRNA standard set of archaea
source("Scripts/standard_set/tRNA_coverage_of_standard_set.R")

# step 3: map out the deviations from the standard set for organisms
source("Scripts/standard_set/deviations_of_standard_set.R")

# step 4: plot the codon coverage of these deviated trna pools to illustrate 
#         impossibilities (codons with no coverage at all)
source("Scripts/standard_set/coverage_plot_deviations.R")

# step 5: compare coverage of this deviation approach with earlier merge attempt
source("Scripts/standard_set/deviation_vs_direct_validation.R")

#-------------------------------------------------------------------------------
# taxonomy and phylogeny
#-------------------------------------------------------------------------------

# Step 1: clean tree labels, identify genomes represented in the phylogeny,
# and create the dataset used for phylogenetic analyses.
source("Scripts/data_prep/taxonomy.R")

# Step 2: create the phylogenetic tree figure and summarize taxonomic
# composition of the genomes retained for phylogenetic analysis.
source("Scripts/creating_plots/tree_plot.R")

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

# step 5: combine the individual amino-acid plots into the final codon-level figure
source("Scripts/creating_plots/combined_codon_analysis_plot.R")

#-------------------------------------------------------------------------------
# codon-pair data preparation
#-------------------------------------------------------------------------------

# Prepare the codon-pair dataset by adding GC3, individual-codon RSCU,
# tRNA coverage states and log2-transformed RDCU.
source("Scripts/data_prep/Codon_pairs.R")

# Run pair-level and constituent-codon analyses
source("Scripts/statistics/codon_pairs_analysis.R")
