review_missing <- function(org_id,
                           missing_df = missing_check,
                           lowq_df = trna_low_quality) {
  
  cat("\n============================\n")
  cat("Organism:", org_id, "\n")
  cat("============================\n\n")
  
  cat("---- Missing check ----\n")
  print(
    missing_df %>%
      filter(organism_id == org_id)
  )
  
  cat("\n---- Low-quality tRNAs ----\n")
  print(
    lowq_df %>%
      filter(organism_id == org_id)
  )
  
  invisible(NULL)
}