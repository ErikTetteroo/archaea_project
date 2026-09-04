# load project
source("R/load_project.R")

# Load data
pairs <- read_csv("Data/cleaned_data/codon_pair_usage_clean.csv")
tree <- read.tree("Data/cleaned_data/cleaned_archaea.tree")

# Keep observations with an assigned coverage state for both codons
pairs_clean <- pairs %>%
  filter(
    !is.na(CM1),
    !is.na(CM2)
  )

# Select codon pairs with sufficient variation in coverage state
pair_state_counts <- pairs_clean %>%
  count(codon_pair, CMP)

pair_state_table <- pair_state_counts %>%
  pivot_wider(
    names_from = CMP,
    values_from = n,
    values_fill = 0
  )

# Retain codon pairs with at least two coverage states containing >=30
# observations.
interesting_pairs <- pair_state_counts %>%
  group_by(codon_pair) %>%
  filter(n >= 30) %>%
  summarise(
    n_states = n(),
    .groups = "drop"
  ) %>%
  filter(n_states >= 2)

interesting_state_table <- pair_state_table %>%
  filter(codon_pair %in% interesting_pairs$codon_pair)

#-------------------------------------------------------------------------------
# Select pairs suitable for the 0/1/2-wobble analysis
#-------------------------------------------------------------------------------

candidate_pairs <- interesting_state_table %>%
  filter(
    GU_X_GU >= 30,
    (GU_X_M + M_X_GU) >= 30,
    M_X_M >= 30
  )

pairs_focus <- pairs_clean %>%
  filter(
    codon_pair %in% candidate_pairs$codon_pair,
    CMP %in% c(
      "M_X_M",
      "GU_X_M",
      "M_X_GU",
      "GU_X_GU"
    )
  ) %>%
  mutate(
    wobbles = case_when(
      CMP == "M_X_M" ~ 0,
      CMP %in% c("GU_X_M", "M_X_GU") ~ 1,
      CMP == "GU_X_GU" ~ 2,
      TRUE ~ NA_real_
    ),
    
    wobble_state = case_when(
      CMP == "M_X_M" ~ "0_GU",
      CMP %in% c("GU_X_M", "M_X_GU") ~ "1_GU",
      CMP == "GU_X_GU" ~ "2_GU",
      TRUE ~ NA_character_
    )
  )

pairs_focus <- pairs_focus %>%
  mutate(
    wobble_state = factor(
      wobble_state,
      levels = c("0_GU", "1_GU", "2_GU")
    )
  )

# Initial descriptive comparison
pair_summary <- pairs_focus %>%
  group_by(codon_pair, wobbles) %>%
  summarise(
    mean_RDCU = mean(log2_RDCU),
    .groups = "drop"
  )

pair_boxplot <- ggplot(
  pair_summary,
  aes(
    x = factor(wobbles),
    y = mean_RDCU
  )
) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x = "Number of GU-wobble-covered codons",
    y = "Mean log2(RDCU)"
  )

pair_boxplot

ggsave(
  filename = "Plots/pair_analysis_wobble_boxplot.png",
  plot = pair_boxplot,
  width = 7,
  height = 5,
  dpi = 300
)

# Phylogenetic models for individual codon pairs
pair_results <- list()

pair_names <- unique(pairs_focus$codon_pair)

for (i in seq_along(pair_names)) {
  
  pair_name <- pair_names[i]
  
  message(
    "Pair ",
    i,
    " / ",
    length(pair_names),
    ": ",
    pair_name
  )
  
  pair_data <- pairs_focus %>%
    filter(codon_pair == pair_name)
  
  rownames(pair_data) <- pair_data$organism
  
  fit_gc <- tryCatch(
    phylolm(
      log2_RDCU ~ GC3,
      phy = tree,
      data = pair_data,
      model = "lambda"
    ),
    error = function(e) {
      message("  GC3 model failed: ", e$message)
      NULL
    }
  )
  
  fit_numeric <- tryCatch(
    phylolm(
      log2_RDCU ~ GC3 + wobbles,
      phy = tree,
      data = pair_data,
      model = "lambda"
    ),
    error = function(e) {
      message("  Numeric wobble model failed: ", e$message)
      NULL
    }
  )
  
  fit_numeric_interaction <- tryCatch(
    phylolm(
      log2_RDCU ~ GC3 + wobbles + GC3:wobbles,
      phy = tree,
      data = pair_data,
      model = "lambda"
    ),
    error = function(e) {
      message("  Numeric interaction model failed: ", e$message)
      NULL
    }
  )
  
  fit_factor <- tryCatch(
    phylolm(
      log2_RDCU ~ GC3 + wobble_state,
      phy = tree,
      data = pair_data,
      model = "lambda"
    ),
    error = function(e) {
      message("  Factor wobble model failed: ", e$message)
      NULL
    }
  )
  
  fit_factor_interaction <- tryCatch(
    phylolm(
      log2_RDCU ~ GC3 + wobble_state + GC3:wobble_state,
      phy = tree,
      data = pair_data,
      model = "lambda"
    ),
    error = function(e) {
      message("  Factor interaction model failed: ", e$message)
      NULL
    }
  )
  
  pair_results[[pair_name]] <- list(
    data = pair_data,
    
    fit_gc = fit_gc,
    fit_num = fit_numeric,
    fit_num_int = fit_numeric_interaction,
    fit_fac = fit_factor,
    fit_fac_int = fit_factor_interaction,
    
    aic = c(
      gc = if (is.null(fit_gc)) NA_real_ else fit_gc$aic,
      numeric = if (is.null(fit_numeric)) NA_real_ else fit_numeric$aic,
      numeric_interaction = if (is.null(fit_numeric_interaction)) {
        NA_real_
      } else {
        fit_numeric_interaction$aic
      },
      factor = if (is.null(fit_factor)) NA_real_ else fit_factor$aic,
      factor_interaction = if (is.null(fit_factor_interaction)) {
        NA_real_
      } else {
        fit_factor_interaction$aic
      }
    )
  )
}

# Summarise pair-level model comparisons
pair_summary_results <- map_dfr(
  names(pair_results),
  function(pair_name) {
    
    aics <- pair_results[[pair_name]]$aic
    
    tibble(
      codon_pair = pair_name,
      gc = unname(aics["gc"]),
      numeric = unname(aics["numeric"]),
      numeric_interaction = unname(aics["numeric_interaction"]),
      factor = unname(aics["factor"]),
      factor_interaction = unname(aics["factor_interaction"])
    )
  }
)

pair_summary_results <- pair_summary_results %>%
  rowwise() %>%
  mutate(
    best = min(
      c(
        gc,
        numeric,
        numeric_interaction,
        factor,
        factor_interaction
      ),
      na.rm = TRUE
    ),
    
    delta_gc = gc - best,
    delta_numeric = numeric - best,
    delta_numeric_interaction = numeric_interaction - best,
    delta_factor = factor - best,
    delta_factor_interaction = factor_interaction - best,
    
    # Improvement over the GC3-only model
    delta_to_gc = gc - min(
      c(
        gc,
        numeric,
        numeric_interaction,
        factor,
        factor_interaction
      ),
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  mutate(
    best_model = case_when(
      delta_gc == 0 ~ "GC3",
      delta_numeric == 0 ~ "GC3 + numeric wobble",
      delta_numeric_interaction == 0 ~
        "GC3 + numeric wobble + interaction",
      delta_factor == 0 ~ "GC3 + factor(wobble)",
      delta_factor_interaction == 0 ~
        "GC3 + factor(wobble) + interaction",
      TRUE ~ NA_character_
    )
  )

# Plot improvement over the GC3-only model
pair_delta_aic_plot <- ggplot(
  pair_summary_results,
  aes(
    x = reorder(codon_pair, delta_to_gc),
    y = delta_to_gc,
    fill = best_model
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Codon pair",
    y = "Improvement in AIC relative to GC3-only model",
    fill = "Best model"
  )

pair_delta_aic_plot

ggsave(
  filename = "Plots/pair_analysis_delta_AIC.png",
  plot = pair_delta_aic_plot,
  width = 8,
  height = 10,
  dpi = 300
)

# Generate plots for the 10 pairs with the greatest improvement
top_pairs <- pair_summary_results %>%
  arrange(desc(delta_to_gc)) %>%
  slice_head(n = 10) %>%
  pull(codon_pair)

dir.create(
  "Plots/pair_plots",
  showWarnings = FALSE,
  recursive = TRUE
)

for (pair_name in top_pairs) {
  
  delta <- pair_summary_results %>%
    filter(codon_pair == pair_name) %>%
    pull(delta_to_gc)
  
  plot <- plot_pair_interaction(
    pair_name,
    pair_results
  ) +
    labs(
      title = paste0(
        pair_name,
        " (ΔAIC = ",
        round(delta, 1),
        ")"
      )
    )
  
  ggsave(
    filename = file.path(
      "Plots/pair_plots",
      paste0(pair_name, ".jpeg")
    ),
    plot = plot,
    width = 6,
    height = 5,
    dpi = 300
  )
}

# Prepare data for constituent-codon follow-up analysis

# Number of synonymous codons per amino acid.
# This is the maximum possible RSCU for a codon if all synonymous codons
# are represented equally in the denominator.
codon_table <- give_codon_table()

aa_sizes <- codon_table %>%
  count(
    aa,
    name = "max_RSCU"
  )

pairs_focus <- pairs_focus %>%
  left_join(
    aa_sizes,
    by = c("aa1" = "aa")
  ) %>%
  rename(max_RSCU1 = max_RSCU) %>%
  left_join(
    aa_sizes,
    by = c("aa2" = "aa")
  ) %>%
  rename(max_RSCU2 = max_RSCU) %>%
  mutate(
    scaled_RSCU1 = RSCU1 / max_RSCU1,
    scaled_RSCU2 = RSCU2 / max_RSCU2
  )

# Avoid logit values of +/- infinity when scaled RSCU is exactly 0 or 1.
eps <- 1e-4

pairs_focus <- pairs_focus %>%
  mutate(
    scaled_RSCU1 = pmin(
      pmax(scaled_RSCU1, eps),
      1 - eps
    ),
    scaled_RSCU2 = pmin(
      pmax(scaled_RSCU2, eps),
      1 - eps
    ),
    
    logit_RSCU1 = qlogis(scaled_RSCU1),
    logit_RSCU2 = qlogis(scaled_RSCU2)
  )

# Constituent-codon follow-up analysis
codon_results <- list()

pair_names <- unique(pairs_focus$codon_pair)

for(i in seq_along(pair_names)) {
  
  pair_name <- pair_names[i]
  
  cat(i, "/", length(pair_names), "-", pair_name, "\n")
  
  pair_dat <- pairs_focus %>%
    filter(codon_pair == pair_name)
  
  ## -----------------------------------------
  ## Determine category
  ## -----------------------------------------
  
  n_M_GU <- sum(pair_dat$CMP == "M_X_GU")
  n_GU_M <- sum(pair_dat$CMP == "GU_X_M")
  
  ## =========================================
  ## Category A
  ## =========================================
  
  if(n_M_GU < 30 & n_GU_M < 30){
    
    cat("  Category A - skipped\n")
    
    next
  }
  
  ## =========================================
  ## Category B
  ## Analyse codon 2
  ## =========================================
  
  if(n_M_GU >= 30 & n_GU_M < 30){
    
    cat("  Category B - codon2\n")
    
    dat <- pair_dat %>%
      filter(CMP != "GU_X_M")
    
    result <- analyze_pair_codon(
      dat = dat,
      tree = tree,
      response = "logit_RSCU2"
    )
    
    result$pair <- pair_name
    result$codon <- 2
    result$response <- "logit_RSCU2"
    result$category <- "B"
    
    codon_results[[paste0(pair_name, "_2")]] <- result
    
    next
  }
  
  ## =========================================
  ## Category C
  ## Analyse codon 1
  ## =========================================
  
  if(n_M_GU < 30 & n_GU_M >= 30){
    
    cat("  Category C - codon1\n")
    
    dat <- pair_dat %>%
      filter(CMP != "M_X_GU")
    
    result <- analyze_pair_codon(
      dat = dat,
      tree = tree,
      response = "logit_RSCU1"
    )
    
    result$pair <- pair_name
    result$codon <- 1
    result$response <- "logit_RSCU1"
    result$category <- "C"
    
    codon_results[[paste0(pair_name, "_1")]] <- result
    
    next
  }
  
  ## =========================================
  ## Category D
  ## Analyse both codons
  ## =========================================
  
  if(n_M_GU >= 30 & n_GU_M >= 30){
    
    cat("  Category D - both codons\n")
    
    ## ---- Codon 1 ----
    
    result1 <- analyze_pair_codon(
      dat = pair_dat,
      tree = tree,
      response = "logit_RSCU1"
    )
    
    result1$pair <- pair_name
    result1$codon <- 1
    result1$response <- "logit_RSCU1"
    result1$category <- "D"
    
    codon_results[[paste0(pair_name, "_1")]] <- result1
    
    
    ## ---- Codon 2 ----
    
    result2 <- analyze_pair_codon(
      dat = pair_dat,
      tree = tree,
      response = "logit_RSCU2"
    )
    
    result2$pair <- pair_name
    result2$codon <- 2
    result2$response <- "logit_RSCU2"
    result2$category <- "D"
    
    codon_results[[paste0(pair_name, "_2")]] <- result2
    
  }
  
}


codon_summary_results <- map_dfr(
  
  names(codon_results),
  
  function(x){
    
    res <- codon_results[[x]]
    
    aics <- res$aic
    
    tibble(
      
      result = x,
      
      pair = res$pair,
      codon = res$codon,
      response = res$response,
      category = res$category,
      
      gc = aics["gc"],
      factor = aics["factor"],
      factor_interaction = aics["factor_interaction"]
      
    )
    
  }
  
)

codon_summary_results <- codon_summary_results %>%
  
  rowwise() %>%
  
  mutate(
    
    best = min(
      c(gc,
        factor,
        factor_interaction),
      na.rm = TRUE
    ),
    
    delta_gc = gc - best,
    delta_factor = factor - best,
    delta_factor_interaction =
      factor_interaction - best
    
  ) %>%
  
  ungroup()

codon_summary_results <- codon_summary_results %>%
  
  mutate(
    
    best_model = case_when(
      
      delta_gc == 0 ~
        "GC3",
      
      delta_factor == 0 ~
        "GC3 + CMP",
      
      delta_factor_interaction == 0 ~
        "GC3 + CMP + GC3:CMP"
      
    )
    
  )

codon_summary_results <- codon_summary_results %>%
  
  mutate(
    
    delta_to_gc =
      gc - pmin(
        gc,
        factor,
        factor_interaction
      )
    
  )

top_codons <- codon_summary_results %>%
  arrange(desc(delta_to_gc)) %>%
  pull(result)

top_codons <- top_codons[1:10]

dir.create(
  "Plots/codon_plots",
  showWarnings = FALSE
)

for(result_name in top_codons){
  
  delta <- codon_summary_results %>%
    filter(result == result_name) %>%
    pull(delta_to_gc)
  
  p <- plot_codon_interaction(
    result_name,
    codon_results
  )
  
  p <- p +
    labs(
      title = paste0(
        result_name,
        " (ΔAIC = ",
        round(delta,1),
        ")"
      )
    )
  
  ggsave(
    file.path(
      "Plots/codon_plots",
      paste0(result_name, ".jpeg")
    ),
    p,
    width = 6,
    height = 5
  )
  
}


codon_heatmap <- codon_summary_results %>%
  select(
    pair,
    codon,
    result,
    delta_to_gc
  ) %>%
  mutate(
    codon_number = case_when(
      grepl("_1$", result) ~ "Codon1",
      grepl("_2$", result) ~ "Codon2"
    )
  ) %>%
  select(pair, codon_number, delta_to_gc)

codon_heatmap <- codon_heatmap %>%
  pivot_wider(
    names_from = codon_number,
    values_from = delta_to_gc
  )

heatmap_dat <- pair_summary_results %>%
  select(
    pair = codon_pair,
    Pair = delta_to_gc
  ) %>%
  left_join(
    codon_heatmap,
    by = "pair"
  )

heatmap_dat <- heatmap_dat %>%
  arrange(desc(Pair))

heatmap_dat$pair <- factor(
  heatmap_dat$pair,
  levels = rev(heatmap_dat$pair)
)

heatmap_long <- heatmap_dat %>%
  pivot_longer(
    cols = c(Pair, Codon1, Codon2),
    names_to = "analysis",
    values_to = "delta_AIC"
  )


constituent_heatmap <- ggplot(
  heatmap_long,
  aes(
    x = analysis,
    y = pair,
    fill = delta_AIC
  )
) +
  geom_tile(
    colour = "white"
  ) +
  geom_text(
    aes(
      label = ifelse(
        is.na(delta_AIC),
        "",
        round(delta_AIC,1)
      )
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    na.value = "grey90"
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL,
    fill = "ΔAIC\nvs GC3"
  )
constituent_heatmap

ggsave(
  filename = "Plots/constituent_codon_heatmap.png",
  plot = constituent_heatmap,
  width = 7,
  height = 10,
  dpi = 300
)
