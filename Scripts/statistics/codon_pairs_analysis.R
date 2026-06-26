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
ggplot(
  pairs_focus,
  aes(factor(wobbles), log2_RDCU)
) +
  geom_boxplot()


# fit simple model
fit <- lmer(
  log2_RDCU ~ GC3 + factor(wobbles) + (1|codon_pair),
  data = pairs_focus
)

summary(fit)

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

test <- pairs_focus[pairs_focus$organism==unique(pairs_focus$organism)[25],]

fit <- lmer(
  log2_RDCU ~ GC3 + factor(wobbles) +
    (1 | codon_pair),
  data = pairs_focus
)

org_resid <- pairs_focus %>%
  mutate(resid = residuals(fit)) %>%
  group_by(organism) %>%
  summarise(
    mean_resid = mean(resid),
    .groups = "drop"
  )

library(phytools)

phylosig(
  tree,
  setNames(
    org_resid$mean_resid,
    org_resid$organism
  ),
  method = "lambda",
  test = TRUE
)

length(organism_summary$organism[!is.na(organism_summary$double_wobble_effect)])

length(pairs_focus$organism[pairs_focus$organism==unique(pairs_focus$organism)[1]])
head(pairs_focus)


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

pair_results$TTGCCG
pair_summary_results[pair_summary_results$codon_pair=="TTGCCG",]

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

length(pair_summary_results$codon_pair[pair_summary_results$delta_factor==0])
length(pair_summary_results$codon_pair[pair_summary_results$delta_numeric==0])
length(pair_summary_results$codon_pair[pair_summary_results$delta_gc==0])

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

pair_summary_results %>%
  filter(delta_factor == 0) %>%
  arrange(delta_gc)

pair_summary_results <- pair_summary_results %>%
  mutate(
    delta_to_gc = gc - pmin(gc, numeric, numeric_interaction ,factor, factor_interaction)
  )

pair_summary_results[pair_summary_results$delta_to_gc>0,c(1,7:9)]

#visualize results

pair_summary_results %>%
  mutate(
    best_model = case_when(
      delta_gc == 0 ~ "GC3 only",
      delta_numeric == 0 ~ "Numeric wobble",
      delta_factor == 0 ~ "Factor wobble"
    )
  ) %>%
  count(best_model)

ggplot(
  pair_summary_results[pair_summary_results$delta_to_gc>0,],
  aes(
    reorder(codon_pair, delta_to_gc),
    delta_to_gc
  )
) +
  geom_col() +
  coord_flip()


ggplot(
  pairs_focus[pairs_focus$codon_pair=="CCGGGG",],
  aes(
    GC3,
    log2_RDCU,
    colour = wobble_state
  )
) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm")

pair_o[1]

pair <- pair_results[["TTGCCG"]]

best_fit <- pair$fit_fac_int   # example
dat <- pair$data

pred_dat <- expand.grid(
  GC3 = seq(
    min(dat$GC3),
    max(dat$GC3),
    length.out = 200
  ),
  wobble_state = levels(dat$wobble_state)
)

X <- model.matrix(
  ~ GC3 + wobble_state + GC3:wobble_state,
  data = pred_dat
)

b <- coef(best_fit)

pred_dat$pred <- as.numeric(
  X %*% b
)

V <- vcov(best_fit)

pred_dat$se <- sqrt(
  diag(
    X %*% V %*% t(X)
  )
)

pred_dat <- pred_dat %>%
  mutate(
    lower = pred - 1.96 * se,
    upper = pred + 1.96 * se
  )

ggplot() +
  
  geom_point(
    data = pair$data,
    aes(
      GC3,
      log2_RDCU,
      color = wobble_state
    ),
    alpha = 0.4
  ) +
  
  geom_ribbon(
    data = pred_dat,
    aes(
      x = GC3,
      ymin = lower,
      ymax = upper,
      fill = wobble_state
    ),
    alpha = 0.2
  ) +
  
  geom_line(
    data = pred_dat,
    aes(
      GC3,
      pred,
      color = wobble_state
    ),
    linewidth = 1
  )

ggplot(pair$data,
       aes(GC3,
           log2_RDCU,
           color = wobble_state)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y ~ x,
              se = FALSE)



plot_pair_interaction <- function(pair_name,
                                  pair_results) {
  
  pair <- pair_results[[pair_name]]
  
  dat <- pair$data
  fit <- pair$fit_fac_int
  
  if(is.null(fit)) {
    stop("No interaction model found for ", pair_name)
  }
  
  pred_dat <- expand.grid(
    GC3 = seq(
      min(dat$GC3),
      max(dat$GC3),
      length.out = 200
    ),
    wobble_state = levels(dat$wobble_state)
  )
  
  X <- model.matrix(
    ~ GC3 * wobble_state,
    data = pred_dat
  )
  
  b <- as.numeric(coef(fit))
  
  pred_dat$pred <- as.numeric(
    X %*% b
  )
  
  V <- vcov(fit)
  
  pred_dat$se <- sqrt(
    diag(
      X %*% V %*% t(X)
    )
  )
  
  pred_dat <- pred_dat %>%
    mutate(
      lower = pred - 1.96 * se,
      upper = pred + 1.96 * se
    )
  
  ggplot() +
    
    geom_point(
      data = dat,
      aes(
        GC3,
        log2_RDCU,
        color = wobble_state
      ),
      alpha = 0.4
    ) +
    
    geom_ribbon(
      data = pred_dat,
      aes(
        x = GC3,
        ymin = lower,
        ymax = upper,
        fill = wobble_state
      ),
      alpha = 0.2,
      colour = NA
    ) +
    
    geom_line(
      data = pred_dat,
      aes(
        GC3,
        pred,
        color = wobble_state
      ),
      linewidth = 1
    ) +
    
    labs(
      title = pair_name,
      y = "log2 RDCU",
      x = "GC3"
    ) +
    
    theme_bw()
}

top_pairs <- pair_summary_results %>%
  arrange(desc(delta_to_gc)) %>%
  slice(1:13) %>%
  pull(codon_pair)

top_pairs <- top_pairs[1:13]

dir.create("pair_plots", showWarnings = FALSE)

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
      "pair_plots",
      paste0(pair_name, ".pdf")
    ),
    p,
    width = 6,
    height = 5
  )
}

library(patchwork)

plots <- lapply(
  top_pairs,
  plot_pair_interaction,
  pair_results = pair_results
)

wrap_plots(plots)


plot_pair_interaction("CCGGGG",pair_results)

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

top_pair_info$codon_pair <- factor(top_pair_info$codon_pair,
                                   levels = top_pairs)

#--------------------------

pair_name <- "ACGGGG"

state_map <- pairs_focus %>%
  filter(codon_pair == pair_name) %>%
  select(
    organism,
    wobble_state
  )

pair_info <- pairs_focus %>%
  filter(codon_pair == pair_name) %>%
  slice(1)

aa1_target <- pair_info$aa1
aa2_target <- pair_info$aa2

syn_pairs <- pairs_c %>%
  filter(
    aa1 == aa1_target,
    aa2 == aa2_target
  )

plot_dat <- syn_pairs %>%
  inner_join(
    state_map,
    by = "organism"
  )

plot_sum <- plot_dat %>%
  group_by(
    wobble_state,
    codon_pair
  ) %>%
  summarise(
    mean_RDCU = mean(log2_RDCU),
    .groups = "drop"
  )

plot_dat <- plot_dat %>%
  mutate(
    focal = codon_pair == pair_name
  )

ggplot(
  plot_dat,
  aes(
    codon_pair,
    log2_RDCU,
    fill = focal
  )
) +
  geom_boxplot() +
  
  facet_wrap(
    ~ wobble_state
  ) +
  
  theme_bw() +
  theme(
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1
      )
  )

ggplot(
  plot_sum,
  aes(
    codon_pair,
    mean_RDCU,
    fill = focal
  )
) +
  geom_boxplot() +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.5
  ) +
  
  facet_wrap(
    ~ wobble_state,
    scales = "free_y"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1
      )
  )
