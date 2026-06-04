local({
  
#load project
source("R/load_project.R")

#create codon_table
codon_table <- give_codon_table()

codons <- codon_table %>%
  distinct(codon_rna) %>%
  pull()

anticodons <- codon_table %>%
  distinct(anticodon) %>%
  pull()

#create empty codon matrix
coverage_mat <- give_empty_codon_matrix(codons = codons,
                        anticodons = anticodons)

#check potential spots for inosine/super wobble
inosine_allowed <- codon_table %>%
  distinct(anticodon) %>%
  mutate(
    inosine_allowed = sapply(
      anticodon,
      can_use_inosine,
      codon_tbl = codon_table
    )
  )

superwobble_allowed <- codon_table %>%
  distinct(anticodon) %>%
  mutate(
    superwobble_allowed = sapply(
      anticodon,
      can_superwobble,
      codon_tbl = codon_table
    )
  )

#fill matrix
for (anticodon in anticodons) {
  
  anti_aa <- codon_table %>%
    filter(anticodon == !!anticodon) %>%
    pull(aa) %>%
    first()
  
  anti_inosine <- inosine_allowed %>%
    filter(anticodon == !!anticodon) %>%
    pull(inosine_allowed)
  
  anti_superwobble <- superwobble_allowed %>%
    filter(anticodon == !!anticodon) %>%
    pull(superwobble_allowed)
  
  for (codon in codons) {
    
    codon_aa <- codon_table %>%
      filter(codon_rna == !!codon) %>%
      pull(aa) %>%
      first()
    
    # only synonymous decoding allowed
    if (anti_aa != codon_aa) {
      next
    }
    
    coverage_mat[anticodon, codon] <- check_pair(
      codon = codon,
      anticodon = anticodon,
      
      inosine_allowed = anti_inosine,
      superwobble_allowed = anti_superwobble,
      
      # organism-specific settings
      allow_inosine = FALSE,     # archaea
      allow_superwobble = TRUE
    )
  }
}

#long format
df_long <- coverage_mat %>%
  as_tibble(rownames = "anticodon") %>%
  pivot_longer(
    cols = -anticodon,
    names_to = "codon",
    values_to = "pairing"
  ) %>%
  mutate(
    pairing = na_if(pairing, "0")
  )

df_long <- df_long %>%
  left_join(
    codon_table %>%
      select(codon_rna, aa),
    by = c("codon" = "codon_rna")
  )

source("R/plot_functions/coverage_plot_alternative.R")

p <- create_alternative_coverage_plot(df_long,
                                 codon_table
                                 )

# save plot
ggsave(
  filename = "Plots/coverage_plot.png",
  plot = p,
  width = 12,
  height = 10,
  dpi = 300
)

# export table
df_out <- df_long %>%
  mutate(
    anticodon = reverse_complement_format(anticodon)
  )

df_out <- df_out %>%
  mutate(
    pairing = if_else(
      anticodon == "CAT" & codon == "AUA" & is.na(pairing),
      "M2",
      pairing
    )
  )

write_csv(
  df_out,
  file = "Data/cleaned_data/coverage_table.csv"
)

})