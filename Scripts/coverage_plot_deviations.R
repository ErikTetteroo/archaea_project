# Archaeal tRNA Coverage Plotting Functions
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)
source("Functions/coverage_plot_function.R")

standard_set_coverage_table <- read_csv("Data/standard_set_coverage_table.csv")
deviating_aa_sets <- read.csv("Data/aa_deviating_state_counts.csv")

# =========================================================
# USER INPUTS
# =========================================================

GUweight <- 0.5
deviation_row <- deviating_aa_sets[1,]

# =========================================================
# CODON ORDER
# =========================================================

codon_order <- list(
  Phe = c("TTT","TTC"),
  Leu = c("TTA","TTG","CTT","CTC","CTA","CTG"),
  Ile = c("ATT","ATC","ATA"),
  Met = c("ATG"),
  Val = c("GTT","GTC","GTA","GTG"),
  Ser = c("TCT","TCC","TCA","TCG","AGT","AGC"),
  Pro = c("CCT","CCC","CCA","CCG"),
  Thr = c("ACT","ACC","ACA","ACG"),
  Ala = c("GCT","GCC","GCA","GCG"),
  Tyr = c("TAT","TAC"),
  His = c("CAT","CAC"),
  Gln = c("CAA","CAG"),
  Asn = c("AAT","AAC"),
  Lys = c("AAA","AAG"),
  Asp = c("GAT","GAC"),
  Glu = c("GAA","GAG"),
  Cys = c("TGT","TGC"),
  Trp = c("TGG"),
  Arg = c("CGT","CGC","CGA","CGG","AGA","AGG"),
  Gly = c("GGT","GGC","GGA","GGG")
)

codon_order_u <- lapply(codon_order, function(x) {
  str_replace_all(x, "T", "U")
})

#plot_deviation_comparison(
#   standard_set_coverage_table,
#   deviation_row,
#   GUweight = GUweight
# )

weighted_standard_table <- standard_set_coverage_table %>%
  mutate(
    GUw = GUw * GUweight,
    CV = M + GUw
  )

# =========================================================
# SAVE ALL DEVIATION PLOTS
# =========================================================

output_base <- "Plots/coverage_plots"

for (i in seq_len(nrow(deviating_aa_sets))) {
  
  deviation_row <- deviating_aa_sets[i, ]
  
  aa <- deviation_row$tRNA_type
  signature <- deviation_row$aa_signature
  freq <- deviation_row$n
  
  # -------------------------------------------------------
  # Create amino-acid-specific output folder
  # -------------------------------------------------------
  
  aa_dir <- file.path(output_base, aa)
  
  if (!dir.exists(aa_dir)) {
    dir.create(aa_dir, recursive = TRUE)
  }
  
  # -------------------------------------------------------
  # Parse codon deviations from signature
  # -------------------------------------------------------
  
  parsed <- str_match_all(
    signature,
    "([AUGC]{3})\\((-?[0-9.]+),(-?[0-9.]+),(-?[0-9.]+)\\)"
  )[[1]]
  
  codons <- parsed[,2]
  delta_cv <- as.numeric(parsed[,5])
  
  # Determine add/del labels
  dev_labels <- ifelse(delta_cv > 0, "add", "del")
  
  # Build codon_dev filename section
  codon_dev_parts <- map2_chr(
    codons,
    dev_labels,
    ~ paste(.x, .y, sep = "_")
  )
  
  codon_dev_string <- paste(codon_dev_parts, collapse = "_")
  
  # -------------------------------------------------------
  # Create filename
  # -------------------------------------------------------
  
  filename <- paste(
    i,
    freq,
    codon_dev_string,
    sep = "_"
  )
  
  filename <- paste0(filename, ".png")
  
  filepath <- file.path(aa_dir, filename)
  
  # -------------------------------------------------------
  # Generate plot
  # -------------------------------------------------------
  
  p <- plot_deviation_comparison(
    weighted_standard_table,
    deviation_row,
    GUweight = GUweight
  )
  
  # -------------------------------------------------------
  # Save plot
  # -------------------------------------------------------
  
  ggsave(
    filename = filepath,
    plot = p,
    width = 10,
    height = 4,
    dpi = 300
  )
  
  message("Saved: ", filepath)
}
