create_alternative_coverage_plot <- function(
    data,
    codon_table,
    coverage_table = NULL,
    aa_subset = NULL,
    title = "coverage of codons by anticodons"
) {
  ###########################################################################
  # INPUT HANDLING
  
  # detect whether input is:
  # 1) full pairing table
  # 2) codon coverage summary table
  
  if (all(c("anticodon", "pairing") %in% colnames(data))) {
    
    # already long format
    df_long <- data
    
    active_anticodons <- unique(df_long$anticodon)
    
  } else if (all(c("codon", "CM") %in% colnames(data))) {
    
    if (is.null(coverage_table)) {
      stop("coverage_table must be provided for summary input.")
    }
    
    #########################################################################
    # reconstruct plotting table from CM categories
    coverage_lookup <- coverage_table #%>%
     # filter(!is.na(pairing))
    
    df_long <- coverage_lookup %>%
      
      left_join(
        data,
        by = "codon"
      ) %>%
      
      group_by(codon, pairing) %>%
      
      mutate(
        keep_rank = row_number(),
        
        keep_n = case_when(
          pairing == "M"  ~ first(M),
          pairing == "M2" ~ first(M2),
          pairing == "GU" ~ first(GUw),
          pairing == "SU" ~ first(SUw),
          TRUE ~ 0
        ),
        
        keep = keep_rank <= keep_n
      ) %>%
      
      ungroup() %>%
      
      mutate(
        pairing = if_else(
          keep,
          pairing,
          NA_character_
        )
      ) %>%
      
      select(-keep_rank, -keep_n, -keep)
    
    df_long$anticodon <- trnascan_to_anticodon(df_long$anticodon)
    
    active_anticodons <- df_long %>%
      filter(!is.na(pairing)) %>%
      pull(anticodon) %>%
      unique()
    
  } else {
    
    stop("Input format not recognized.")
    
  }
  
  ###########################################################################
  # OPTIONAL AA SUBSETTING
  
  if (!is.null(aa_subset)) {
    
    valid_codons <- codon_table %>%
      filter(aa == aa_subset) %>%
      pull(codon_rna)
    
    valid_anticodons <- codon_table %>%
      filter(aa == aa_subset) %>%
      pull(anticodon)
    
    df_long <- df_long %>%
      filter(
        codon %in% valid_codons,
        anticodon %in% valid_anticodons
      )
  }
  
  ###########################################################################
  # ORDERING
  
  codons <- codon_table %>%
    distinct(codon_rna) %>%
    pull(codon_rna)
  
  # anticodons ordered as complements to codons
  anticodons <- comp(codons)
  
  df_long <- df_long %>%
    
#    mutate(
#      anticodon_present = anticodon %in% active_anticodons
#    ) %>%
    
    mutate(
      codon = factor(
        codon,
        levels = rev(codons)
      ),
      
      anticodon = factor(
        anticodon,
        levels = anticodons
      )
    )
  
  ###########################################################################
  # AA BOUNDARIES
  
  aa_boundaries <- codon_table %>%
    distinct(codon_rna, aa) %>%
    mutate(
      codon_rna = factor(codon_rna, levels = rev(codons))
    ) %>%
    arrange(codon_rna) %>%
    mutate(position = row_number()) %>%
    group_by(aa) %>%
    summarise(
      start = min(position),
      end = max(position),
      middle = mean(c(start, end)),
      .groups = "drop"
    )
  
  ###########################################################################
  # PLOTTING
  
  ggplot(
    df_long,
    aes(
      x = anticodon,
      y = codon
    )
  ) +
    
    # grey background for absent anticodons
    geom_tile(color = "grey90") +
    
    # actual pairings
    geom_tile(
      aes(fill = pairing),
      color = "grey90"
    ) +
    
    geom_hline(
      data = aa_boundaries,
      aes(yintercept = end + 0.5),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    
    geom_text(
      data = aa_boundaries,
      aes(
        x = 1,
        y = middle,
        label = aa
      ),
      inherit.aes = FALSE,
      hjust = 1,
      size = 3
    ) +
    
    scale_fill_manual(
      values = c(
        "M"  = "#2166AC",
        "M2" = "#542788",
        "GU" = "#67A9CF",
        "I"  = "#FDAE61",
        "SU" = "#B2182B"
      ),
      na.value = "white"
    ) +
    
    coord_cartesian(clip = "off") +
    
    theme_minimal() +
    
    theme(
      panel.grid = element_blank(),
      
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        size = 6,
        
#        color = ifelse(
#          anticodons %in% active_anticodons,
#          "black",
#          "grey70"
#        )
      ),
      
      axis.text.y = element_text(size = 6),
      
      plot.margin = margin(
        t = 10,
        r = 10,
        b = 10,
        l = 40
      )
    ) +
    
    labs(
      fill = "Pairing type",
      x = "Anticodon",
      y = "Codon"
    ) +
    
    ggtitle(title)
}

create_alternative_coverage_plot(standard_set_coverage_table,
                                 codon_table,
                                 coverage_table)

