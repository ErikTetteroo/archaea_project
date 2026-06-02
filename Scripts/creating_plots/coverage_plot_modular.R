source("R/load_project.R")

standard_set_coverage_table <- read_csv("Data/standard_set_coverage_table.csv")

# =========================================================
# USER OPTIONS
# =========================================================

GUweight <- 0.5

# Set to NULL for all amino acids
# Or specify one amino acid, e.g. "Arg"
selected_aa <- 'Leu'

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
  Gly = c("GGT","GGC","GGA","GGG"),
  TER = c("TAA","TAG","TGA")
)

# =========================================================
# PREP CODON ORDER
# =========================================================

# Convert T -> U
codon_order_u <- lapply(codon_order, function(x) {
  str_replace_all(x, "T", "U")
})

# Remove TER from AA grouping
codon_order_u$TER <- NULL

# Optionally subset to one amino acid
if (!is.null(selected_aa)) {
  codon_order_u <- codon_order_u[selected_aa]
}

# Flatten codon vector
codon_levels <- unlist(codon_order_u)

# Add UGA only for full plot
if (is.null(selected_aa)) {
  codon_levels <- c(codon_levels, "UGA")
}

# =========================================================
# AMINO ACID LABELS / BOUNDARIES
# =========================================================

aa_df <- tibble(
  amino_acid = names(codon_order_u),
  codons = codon_order_u
) %>%
  unnest(codons) %>%
  mutate(x = row_number())

aa_labels <- aa_df %>%
  group_by(amino_acid) %>%
  summarise(
    xmin = min(x),
    xmax = max(x),
    xmid = mean(c(xmin, xmax)),
    .groups = "drop"
  )

aa_boundaries <- aa_labels %>%
  mutate(boundary = xmax + 0.5) %>%
  pull(boundary)

# =========================================================
# PREP COVERAGE TABLE
# =========================================================

standard_set_coverage_table <- standard_set_coverage_table %>%
  mutate(
    GUw = GUweight * GUw
  )

# Add empty UGA row only for full plot
if (is.null(selected_aa)) {
  
  plot_input <- standard_set_coverage_table %>%
    bind_rows(
      tibble(
        codon = "UGA",
        ile2_AUA = FALSE,
        M = 0,
        GUw = 0,
        CV = 0
      )
    ) %>%
    distinct(codon, .keep_all = TRUE)
  
} else {
  
  plot_input <- standard_set_coverage_table
}

# Keep only relevant codons
plot_input <- plot_input %>%
  filter(codon %in% codon_levels)

# =========================================================
# RESHAPE FOR PLOTTING
# =========================================================

plot_df <- plot_input %>%
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

# =========================================================
# PLOT
# =========================================================

ggplot(plot_df, aes(x = codon, y = coverage, fill = fill_group)) +
  
  geom_col(position = position_stack(reverse = TRUE)) +
  
  # AA separator lines only when plotting all amino acids
  {
    if (is.null(selected_aa))
      geom_vline(
        xintercept = aa_boundaries,
        color = "grey70",
        linewidth = 0.4
      )
  } +
  
  # AA labels
  geom_text(
    data = aa_labels,
    aes(x = xmid, y = 1.7, label = amino_acid),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  
  labs(
    title = ifelse(
      is.null(selected_aa),
      "Codon Coverage in the Standard Archaeal tRNA Set",
      paste0(selected_aa, " Codon Coverage")
    ),
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
