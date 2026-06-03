create_alternative_coverage_plot <- function(
    df_long,
    codon_table,
    aa = NULL,
    title = "coverage of codons by anticodons"
) {
  
  ###########################################################################
  # optional amino acid subsetting
  
  if (!is.null(aa)) {
    
    valid_codons <- codon_table %>%
      filter(aa == !!aa) %>%
      pull(codon_rna)
    
    valid_anticodons <- codon_table %>%
      filter(aa == !!aa) %>%
      pull(anticodon)
    
    df_long <- df_long %>%
      filter(
        codon %in% valid_codons,
        anticodon %in% valid_anticodons
      )
  }
  
  ###########################################################################
  # ordering
  
  codons <- df_long %>%
    distinct(codon) %>%
    pull(codon)
  
  anticodons <- df_long %>%
    distinct(anticodon) %>%
    pull(anticodon)
  
  df_long <- df_long %>%
    mutate(
      codon = factor(codon, levels = rev(codons)),
      anticodon = factor(anticodon, levels = anticodons)
    )
  
  ###########################################################################
  # amino acid boundaries
  
  aa_boundaries <- codon_table %>%
    filter(codon_rna %in% codons) %>%
    distinct(codon_rna, aa) %>%
    mutate(
      codon_rna = factor(codon_rna, levels = rev(codons))
    ) %>%
    arrange(codon_rna) %>%
    mutate(
      position = row_number()
    ) %>%
    group_by(aa) %>%
    summarise(
      start = min(position),
      end = max(position),
      middle = mean(c(start, end)),
      .groups = "drop"
    )
  
  ###########################################################################
  # plotting
  
  ggplot(
    df_long,
    aes(
      x = anticodon,
      y = codon,
      fill = pairing
    )
  ) +
    
    geom_tile(color = "grey90") +
    
    geom_hline(
      data = aa_boundaries,
      aes(yintercept = end + 0.5),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    
    geom_text(
      data = aa_boundaries,
      aes(
        x = 4,
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
        size = 6
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