#questionable_deviations <- aa_delta_state_counts[aa_delta_state_counts$has_manual==0,]
#questionable_deviations$impossible <- 0

#rown <- 44
#questionable_deviations[rown,]
#questionable_deviations[rown,5] <- 1

#deviations_observerd_by_manuscript <- aa_delta_state_counts[aa_delta_state_counts$has_manual==1,]
#deviations_observerd_by_manuscript$impossible <- 0
#rown <- 33
#deviations_observerd_by_manuscript[rown,]
#deviations_observerd_by_manuscript[rown,5] <- 1
questionable_deviations_u <- read_csv("Data/questionable deviations.csv")
questionable_deviations_u_m <- read_csv("Data/questionable deviations_m.csv")


aa_delta_state_counts <- aa_delta_states %>%
  group_by(tRNA_type, aa_signature) %>%
  summarise(
    n = n(),
    
    # keep all organism IDs as a list-column
    organism_ids = list(unique(organism_id)),
    
    .groups = "drop"
  ) %>%
  arrange(desc(n))

possible_aa_delta_states <- aa_delta_state_counts %>%
  left_join(
    bind_rows(
      questionable_deviations_u,
      questionable_deviations_u_m
    ) %>%
      select(tRNA_type, aa_signature, impossible) %>%
      distinct(),
    
    by = c("tRNA_type", "aa_signature")
  ) 