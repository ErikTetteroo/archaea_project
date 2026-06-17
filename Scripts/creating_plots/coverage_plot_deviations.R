# Archaeal tRNA Coverage Plotting Functions
source("R/load_project.R")

standard_set_coverage_table <- read_csv("Data/manuscript/standard_set_coverage_table.csv")
deviating_aa_sets <- read_csv("Data/manuscript/aa_deviating_state_counts.csv")
coverage_table <- read_csv("Data/cleaned_data/coverage_table.csv")
codon_table <- give_codon_table()


# =========================================================
# SAVE ALL DEVIATION PLOTS
# =========================================================

#remove Sec/Pyl deviations
deviating_aa_sets <- deviating_aa_sets %>%
  filter(!tRNA_type == "SeC",
         !tRNA_type == "Pyl",
         !tRNA_type == "Sup")


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
    "([AUGC]{3})\\((-?[0-9.]+),(-?[0-9.]+),(-?[0-9.]+),(-?[0-9.]+)\\)"
  )[[1]]
  
  codons <- parsed[,2]
  delta_cm <- c(as.numeric(parsed[,3]),
                           as.numeric(parsed[,4]),
                           as.numeric(parsed[,5]),
                           as.numeric(parsed[,6]))
  
  # Determine add/del labels
  dev_labels <- ifelse(sum(delta_cm) > 0, "add", "del")
  
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
  
  
  deviated_coverage_table <- apply_deviation(standard_set_coverage_table,deviation_row$aa_signature)
  
  p1 <- create_alternative_coverage_plot(
    standard_set_coverage_table,
    codon_table,
    coverage_table,
    aa_subset = aa,
    title = "standard coverage"
  )
  
  p2 <- create_alternative_coverage_plot(
    deviated_coverage_table,
    codon_table,
    coverage_table,
    aa_subset = aa,
    title = "deviated coverage"
  )
  
  p <- p1 + p2
  
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

