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

source("R/plot_functions/alternative_coverage_plot.R")