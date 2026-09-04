#load project
source("R/load_project.R")

# Individual amino-acid plots are created in the fold-specific statistics
# scripts and stored in the global environment for use here.

codon_plots <- list(
  Lys = plys,
  Gln = pgln,
  Glu = pglu,
  Val = pval,
  Thr = pthr,
  Ala = pala,
  Gly = pgly,
  Ser = pser,
  Pro = ppro,
  Leu = pleu,
  Arg = parg
)

# Shared colour scale for tRNA coverage states
state_colors <- scale_color_manual(
  name = "State",
  values = c(
    M     = "#1b9e77",
    GU    = "#d95f02",
    GU_SU = "#7570b3",
    GU_GU = "#e7298a"
  ),
  breaks = c("M", "GU", "GU_SU", "GU_GU"),
  drop = TRUE
)

# Apply common formatting to all plots
codon_plots_formatted <- lapply(
  codon_plots,
  function(p) {
    p +
      coord_cartesian(ylim = c(-2.5, 2.5)) +
      state_colors
  }
)

# Combine into a single figure
combined_codon_plot <- wrap_plots(
  codon_plots_formatted,
  ncol = 3
) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

# Display figure
combined_codon_plot

# Save figure
ggsave(
  filename = "Plots/combined_codon_analysis_plot.png",
  plot = combined_codon_plot,
  width = 12,
  height = 10,
  dpi = 300
)