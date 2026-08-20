#load project
source("R/load_project.R")

#read data
pairs <- read_csv("Data/Raw_data/codon_pair_usage_clean.csv")
tree <- read.tree("Data/cleaned_data/cleaned_archaea.tree")

#filter out na's
pairs_c <- pairs %>%
  filter(
    !is.na(CM1),
    !is.na(CM2)
  )

# visualize pairs that have sufficient varying states
pair_states <- pairs_c %>%
  count(codon_pair, CMP) %>%
  arrange(codon_pair, desc(n))

pair_state_table <- pairs_c %>%
  count(codon_pair, CMP) %>%
  tidyr::pivot_wider(
    names_from = CMP,
    values_from = n,
    values_fill = 0
  )

# filter 
interesting_pairs <- pairs_c %>%
  count(codon_pair, CMP) %>%
  group_by(codon_pair) %>%
  filter(n >= 30) %>%          
  summarise(
    states = n(),
    .groups = "drop"
  ) %>%
  filter(states >= 2)

int_pair_state_table <- pair_state_table[pair_state_table$codon_pair %in% interesting_pairs$codon_pair,]


# synergistic wobble effect analysis
candidate_pairs <- int_pair_state_table %>%
  filter(
    GU_X_GU >= 30,
    (GU_X_M + M_X_GU) >= 30,
    M_X_M >= 30
  )

pairs_focus <- pairs_c %>%
  filter(
    codon_pair %in% candidate_pairs$codon_pair,
    CMP %in% c(
      "M_X_M",
      "GU_X_M",
      "M_X_GU",
      "GU_X_GU"
    )
  )

pairs_focus <- pairs_focus %>%
  mutate(
    wobbles =
      (CM1 == "GU") +
      (CM2 == "GU")
  )


# visualize
pairs_focus %>%
  group_by(wobbles) %>%
  summarise(
    mean = mean(log2_RDCU),
    median = median(log2_RDCU),
    sd = sd(log2_RDCU)
  )

pair_summary <- pairs_focus %>%
  group_by(codon_pair, wobbles) %>%
  summarise(
    mean_RDCU = mean(log2_RDCU),
    .groups = "drop"
  )

ggplot(
  pair_summary,
  aes(factor(wobbles), mean_RDCU)
) +
  geom_boxplot()


pair_effects <- pairs_focus %>%
  group_by(codon_pair, CMP) %>%
  summarise(
    mean_RDCU = mean(log2_RDCU),
    n = n(),
    .groups = "drop"
  )

#simple model
fit <- lmer(
  log2_RDCU ~ GC3 + factor(wobbles) +
    (1 | codon_pair) + (1 | organism),
  data = pairs_focus
)

summary(fit)

#factor and numeric version for wobbles 
pairs_focus <- pairs_focus %>%
  mutate(
    wobble_state = case_when(
      CMP == "M_X_M" ~ "0_GU",
      CMP %in% c("GU_X_M", "M_X_GU") ~ "1_GU",
      CMP == "GU_X_GU" ~ "2_GU"
    ),
    
    wobbles = case_when(
      CMP == "M_X_M" ~ 0,
      CMP %in% c("GU_X_M", "M_X_GU") ~ 1,
      CMP == "GU_X_GU" ~ 2
    )
  )




# run individual model for each pair
pair_results <- list()

pair_names <- unique(pairs_focus$codon_pair)

for(i in seq_along(pair_names)) {
  
  pair_name <- pair_names[i]
  
  cat(i, "/", length(pair_names), "-", pair_name, "\n")
  
  pair_dat <- pairs_focus %>%
    filter(codon_pair == pair_name)
  
  
  pair_dat$wobble_state <- factor(
    pair_dat$wobble_state,
    levels = c("0_GU","1_GU","2_GU")
  )
  
  row.names(pair_dat) <- pair_dat$organism
  
  fit_gc <- tryCatch(
    
    phylolm(
      log2_RDCU ~ GC3,
      phy = tree,
      data = pair_dat,
      model = "lambda"
    ),
    
    error = function(e) {
      message("GC model failed for ", pair_name)
      message(e$message)
      NULL
    }
  )
  
  fit_num <- tryCatch(
    
    phylolm(
      log2_RDCU ~ GC3 + wobbles,
      phy = tree,
      data = pair_dat,
      model = "lambda"
    ),
    
    error = function(e) {
      message("State model failed for ", pair_name)
      message(e$message)
      NULL
    }
  )
  
  fit_num_int <- tryCatch(
    
    phylolm(
      log2_RDCU ~ GC3 + wobbles + GC3:wobbles,
      phy = tree,
      data = pair_dat,
      model = "lambda"
    ),
    
    error = function(e) {
      message("State model failed for ", pair_name)
      message(e$message)
      NULL
    }
  )
  
  fit_fac <- tryCatch(
    
    phylolm(
      log2_RDCU ~ GC3 + wobble_state,
      phy = tree,
      data = pair_dat,
      model = "lambda"
    ),
    
    error = function(e) {
      message("State model failed for ", pair_name)
      message(e$message)
      NULL
    }
  )
  
  fit_fac_int <- tryCatch(
    
    phylolm(
      log2_RDCU ~ GC3 + wobble_state + GC3:wobble_state,
      phy = tree,
      data = pair_dat,
      model = "lambda"
    ),
    
    error = function(e) {
      message("State model failed for ", pair_name)
      message(e$message)
      NULL
    }
  )
  
  pair_results[[pair_name]] <- list(
    data = pair_dat,
    
    fit_gc = fit_gc,
    fit_num = fit_num,
    fit_num_int = fit_num_int,
    fit_fac = fit_fac,
    fit_fac_int = fit_fac_int,
    
    aic = c(
      gc = ifelse(is.null(fit_gc), NA, fit_gc$aic),
      numeric = ifelse(is.null(fit_num), NA, fit_num$aic),
      numeric_interaction = ifelse(is.null(fit_num_int), NA, fit_num_int$aic),
      factor = ifelse(is.null(fit_fac), NA, fit_fac$aic),
      factor_interaction = ifelse(is.null(fit_fac_int), NA, fit_fac_int$aic)
    )
  )
}

pair_summary_results <- map_dfr(
  names(pair_results),
  function(x){
    
    aics <- pair_results[[x]]$aic
    
    tibble(
      codon_pair = x,
      gc = aics["gc"],
      numeric = aics["numeric"],
      numeric_interaction = aics["numeric_interaction"],
      factor = aics["factor"],
      factor_interaction = aics["factor_interaction"],
    )
  }
)

pair_summary_results <- pair_summary_results %>%
  rowwise() %>%
  mutate(
    best = min(c(gc, numeric, numeric_interaction, factor, factor_interaction), na.rm = TRUE),
    delta_gc = gc - best,
    delta_numeric = numeric - best,
    delta_numeric_interaction = numeric_interaction - best,
    delta_factor = factor - best,
    delta_factor_interaction = factor_interaction - best
  )

pair_summary_results %>%
  mutate(
    best_model = case_when(
      delta_gc == 0 ~ "GC",
      delta_numeric == 0 ~ "numeric",
      delta_factor == 0 ~ "factor",
      delta_numeric_interaction == 0 ~ "numeric_interaction",
      delta_factor_interaction == 0 ~ "factor_interaction"
    )
  ) %>%
  count(best_model)

pair_summary_results <- pair_summary_results %>%
  mutate(
    delta_to_gc = gc - pmin(gc, numeric, numeric_interaction ,factor, factor_interaction)
  )

#visualize results
pair_summary_results <- pair_summary_results %>%
  mutate(
    best_model = case_when(
      delta_gc == 0 ~ "GC3",
      delta_numeric == 0 ~ "GC3 + num(wobble)",
      delta_numeric_interaction == 0 ~ "GC3 + num(wobble) + GC3:num(wobble)",
      delta_factor == 0 ~ "GC3 + factor(wobble)",
      delta_factor_interaction == 0 ~ "GC3 + factor(wobble) + GC3:factor(wobble)"
    )
  )


ggplot(
  pair_summary_results,
  aes(
    reorder(codon_pair, delta_to_gc),
    delta_to_gc,
    fill = best_model
  )
) +
  geom_col() +
  coord_flip()



# create plots and analyse individual pairs
top_pairs <- pair_summary_results %>%
  arrange(desc(delta_to_gc)) %>%
  pull(codon_pair)

top_pairs <- top_pairs[1:10]

dir.create("Plots/pair_plots", showWarnings = FALSE)

for(pair_name in top_pairs) {
  
  delta <- pair_summary_results %>%
    filter(codon_pair == pair_name) %>%
    pull(delta_to_gc)
  
  p <- plot_pair_interaction(
    pair_name,
    pair_results
  )
  
  p <- p +
    labs(
      title = paste0(
        pair_name,
        "  (ΔAIC = ",
        round(delta,1),
        ")"
      )
    )
  
  ggsave(
    file.path(
      "Plots/pair_plots",
      paste0(pair_name, ".jpeg")
    ),
    p,
    width = 6,
    height = 5
  )
}

top_pair_info <- pairs_focus %>%
  filter(codon_pair %in% top_pairs) %>%
  group_by(codon_pair) %>%
  summarise(
    c1 = first(c1),
    c2 = first(c2),
    aa1 = first(aa1),
    aa2 = first(aa2),
    .groups = "drop"
  )

# visualize distribution of samples within 1 wobble state
pair_state_counts <-
  pairs_focus %>%
  count(codon_pair, CMP) %>%
  tidyr::pivot_wider(
    names_from = CMP,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    one_total = GU_X_M + M_X_GU
  )


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

ggplot(
  codon_summary_results,
  aes(
    reorder(result, delta_to_gc),
    delta_to_gc,
    fill = best_model
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~category, scales = "free_y")

ggplot(
  codon_summary_results,
  aes(
    reorder(result, delta_to_gc),
    delta_to_gc,
    fill = best_model
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~codon)


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

table <- give_codon_table()

# Number of synonymous codons per amino acid
aa_sizes <- table %>%
  count(aa, name = "max_RSCU")

# Join for codon 1
pairs_focus <- pairs_focus %>%
  left_join(
    aa_sizes,
    by = c("aa1" = "aa")
  ) %>%
  rename(max_RSCU1 = max_RSCU)

# Join for codon 2
pairs_focus <- pairs_focus %>%
  left_join(
    aa_sizes,
    by = c("aa2" = "aa")
  ) %>%
  rename(max_RSCU2 = max_RSCU)

# Scale to 0-1
pairs_focus <- pairs_focus %>%
  mutate(
    scaled_RSCU1 = RSCU1 / max_RSCU1,
    scaled_RSCU2 = RSCU2 / max_RSCU2
  )

summary(pairs_focus$scaled_RSCU1)
summary(pairs_focus$scaled_RSCU2)

eps <- 1e-4

pairs_focus <- pairs_focus %>%
  mutate(
    scaled_RSCU1 = pmin(pmax(scaled_RSCU1, eps), 1 - eps),
    scaled_RSCU2 = pmin(pmax(scaled_RSCU2, eps), 1 - eps),
    
    logit_RSCU1 = qlogis(scaled_RSCU1),
    logit_RSCU2 = qlogis(scaled_RSCU2)
  )

ggplot(pairs_focus,
       aes(GC3,
           fill = wobble_state)) +
  geom_histogram(binwidth=.02,
                 position="fill")


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


ggplot(
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


ccgaag_dat <- pairs_focus %>%
  filter(codon_pair == "CCGAAG")

tip_groups <- setNames(
  ccgaag_dat$CMP,
  ccgaag_dat$organism
)

tip_state <- tip_groups[tree$tip.label]

cols <- c(
  "GU_X_GU" = "#d73027",
  "GU_X_M"  = "#fc8d59",
  "M_X_GU"  = "#91bfdb",
  "M_X_M"   = "#4575b4"
)

tip_col <- cols[tip_state]
tip_col[is.na(tip_col)] <- "grey80"

plot(
  tree,
  tip.color = tip_col,
  cex = 0.4
)

legend(
  "topleft",
  legend = names(cols),
  col = cols,
  pch = 19,
  bty = "n"
)


