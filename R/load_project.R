#------------------------------
# Load Packages
#------------------------------

library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(purrr)
library(ape)
library(phylolm)

#------------------------------
# Load basic functions
#------------------------------

source("R/basic/dna_manipulations.R")
source("R/basic/codon_constants.R")


#------------------------------
# Load coverage
#------------------------------

source("R/coverage/matching.R")

#------------------------------
# Load plot functions
#------------------------------

source("R/plot_functions/coverage_plot_alternative.R")

#------------------------------
# Load statistics functions
#------------------------------

source("R/statistics/statistics_functions.R")