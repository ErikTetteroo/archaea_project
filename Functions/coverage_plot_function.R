make_plot_df <- function(coverage_table, selected_aa, GUweight = 0.5) {
  
  codon_levels <- unlist(codon_order_u[selected_aa])
  
  coverage_table <- coverage_table%>%
    filter(codon %in% codon_levels)
  
  plot_df <- coverage_table %>%
    select(codon, ile2_AUA, M, GUw) %>%
    pivot_longer(
      cols = c(M, GUw),
      names_to = "pair_type",
      values_to = "coverage"
    ) %>%
    mutate(
      codon = factor(codon, levels = codon_levels),
      
      fill_group = case_when(
        codon == "AUA" & ile2_AUA ~ "Ile2",
        pair_type == "M" ~ "Exact",
        pair_type == "GUw" ~ "Wobble"
      ),
      
      fill_group = factor(
        fill_group,
        levels = c("Exact", "Ile2", "Wobble")
      )
    )
  
  return(plot_df)
}

make_coverage_plot <- function(plot_df, title_text) {
  
  ggplot(plot_df, aes(x = codon, y = coverage, fill = fill_group)) +
    
    geom_col(position = position_stack(reverse = TRUE)) +
    
    labs(
      title = title_text,
      x = "Codon",
      y = "Coverage Value",
      fill = "Coverage Type"
    ) +
    
    scale_fill_manual(
      values = c(
        "Exact" = "steelblue",
        "Wobble" = "orange",
        "Ile2" = "firebrick"
      ),
      labels = c(
        "Exact" = "Exact match",
        "Wobble" = "GU wobble",
        "Ile2" = "Ile2-mediated AUA decoding"
      )
    ) +
    
    coord_cartesian(
      ylim = c(0, 1.8),
      clip = "off"
    ) +
    
    theme_bw() +
    
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      
      plot.margin = margin(10, 10, 30, 10)
    )
}

# =========================================================
# FUNCTION: APPLY DEVIATION
# =========================================================

apply_deviation <- function(coverage_table, deviation_row) {
  
  signature <- deviation_row$aa_signature
  
  # Parse entries like GCG(-1,0,-1)
  parsed <- str_match_all(
    signature,
    "([AUGC]{3})\\((-?[0-9.]+),(-?[0-9.]+),(-?[0-9.]+)\\)"
  )[[1]]
  
  delta_df <- tibble(
    codon = parsed[,2],
    delta_M = as.numeric(parsed[,3]),
    delta_GUw = as.numeric(parsed[,4]),
    delta_CV = as.numeric(parsed[,5])
  )
  
  updated_table <- coverage_table %>%
    left_join(delta_df, by = "codon") %>%
    mutate(
      delta_M = replace_na(delta_M, 0),
      delta_GUw = replace_na(delta_GUw, 0),
      
      M = pmax(M + delta_M, 0),
      GUw = pmax(GUw + delta_GUw, 0),
      
      CV = M + GUw
    ) %>%
    select(-delta_M, -delta_GUw, -delta_CV)
  
  return(updated_table)
}

# =========================================================
# FUNCTION: COMPARE STANDARD VS DEVIATING STATE
# =========================================================

plot_deviation_comparison <- function(
    standard_table,
    deviation_row,
    GUweight = 0.5
) {
  
  selected_aa <- deviation_row$tRNA_type
  
  # Standard plot dataframe
  standard_plot_df <- make_plot_df(
    standard_table,
    selected_aa,
    GUweight
  )
  
  # Apply deviation
  deviated_table <- apply_deviation(
    standard_table,
    deviation_row
  )
  
  # Deviated plot dataframe
  deviated_plot_df <- make_plot_df(
    deviated_table,
    selected_aa,
    GUweight
  )
  
  p1 <- make_coverage_plot(
    standard_plot_df,
    paste0(selected_aa, " Standard Coverage")
  )
  
  p2 <- make_coverage_plot(
    deviated_plot_df,
    paste0(selected_aa, " Deviated Coverage")
  )
  
  p1 + p2
}