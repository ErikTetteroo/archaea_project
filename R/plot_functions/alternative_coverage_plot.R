create_alternative_coverage_plot <- function(
    df_long,
    title = "coverage of codons by anticodons"
) {
  df_long <- df_long %>%
    mutate(
      codon = factor(codon, levels = rev(codons)),
      anticodon = factor(anticodon, levels = anticodons)
    )
  
  aa_boundaries <- codon_table %>%
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
  
  ggplot(
    df_long,
    aes(
      x = anticodon,
      y = codon,
      fill = pairing
    )
  ) +
    
    geom_tile(color = "grey90") +
    
    # amino acid boundary lines
    geom_hline(
      data = aa_boundaries,
      aes(yintercept = end + 0.5),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    
    # amino acid labels
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