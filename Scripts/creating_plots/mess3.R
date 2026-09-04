library(dplyr)
library(ggplot2)
library(patchwork)

plot_pair_interaction <- function(pair_name,
                                  pair_results) {
  
  pair <- pair_results[[pair_name]]
  
  dat <- pair$data
  fit <- pair$fit_fac_int
  
  if (is.null(fit)) {
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
  
  pred_dat$pred <- as.numeric(X %*% b)
  
  V <- vcov(fit)
  
  pred_dat$se <- sqrt(
    diag(X %*% V %*% t(X))
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
      linewidth = 0.8
    ) +
    labs(
      title = pair_name,
      x = "GC3",
      y = "log2 RDCU"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7),
      legend.position = "none"
    )
}

top_plots <- lapply(seq_len(nrow(top_pairs)), function(i) {
  
  pair_name <- top_pairs$codon_pair[i]
  delta <- top_pairs$delta_to_gc[i]
  
  plot_pair_interaction(
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
})

combined_plot <- wrap_plots(
  top_plots,
  ncol = 2,
  guides = "collect"
) &
  theme(
    legend.position = "bottom"
  )

combined_plot

ggsave(
  "Plots/top_10_pair_interactions.png",
  combined_plot,
  width = 12,
  height = 20,
  dpi = 300
)

ptest <- pair_summary_results[,c(1,13)]

# Get the 10 pairs with the largest AIC improvement
top_pairs <- pair_summary_results %>%
  arrange(desc(delta_to_gc)) %>%
  slice_head(n = 10)

top_pairs <- top_pairs[1:10,]

# Create the plots and store them in a list
top_plots <- vector("list", nrow(top_pairs))

for (i in seq_len(nrow(top_pairs))) {
  
  pair_name <- top_pairs$codon_pair[i]
  delta <- top_pairs$delta_to_gc[i]
  
  p <- plot_pair_interaction(
    pair_name,
    pair_results
  ) +
    labs(
      title = paste0(
        pair_name,
        "  (ΔAIC = ",
        round(delta, 1),
        ")"
      )
    )
  
  top_plots[[i]] <- p
}

# Combine into a 2 x 5 grid
combined_plot <- wrap_plots(
  top_plots,
  ncol = 2
)

combined_plot





library(dplyr)

model_map <- c(
  "GC3" = "fit_gc",
  "GC3 + num(wobble)" = "fit_num",
  "GC3 + num(wobble) + GC3:num(wobble)" = "fit_num_int",
  "GC3 + factor(wobble)" = "fit_fac",
  "GC3 + factor(wobble) + GC3:factor(wobble)" = "fit_fac_int"
)

pair_appendix <- pair_summary_results %>%
  select(
    codon_pair,
    best_model,
    delta_to_gc
  ) %>%
  rowwise() %>%
  mutate(
    model_object = model_map[best_model],
    
    R2 = if (
      is.na(model_object) ||
      is.null(pair_results[[codon_pair]][[model_object]])
    ) {
      NA_real_
    } else {
      pair_results[[codon_pair]][[model_object]]$adj.r.squared
    }
  ) %>%
  ungroup() %>%
  select(
    Codon_pair = codon_pair,
    Best_model = best_model,
    Delta_AIC_to_GC3 = delta_to_gc,
    R2 = R2
  ) %>%
  mutate(
    Delta_AIC_to_GC3 = round(Delta_AIC_to_GC3, 2),
    R2 = round(R2, 3)
  ) %>%
  arrange(desc(Delta_AIC_to_GC3))

pair_appendix
pair_appendix

pair_results$CCGGGG$fit_fac_int$r.squared


head(codon_summary_results)
head(codon_summary_results[,10:14])

codon_results$TTGCCG_2$

fit_fac_int
fit_fac
fit_gc

library(dplyr)

codon_appendix <- codon_summary_results %>%
  select(
    result,
    pair,
    codon,
    best_model,
    delta_to_gc
  ) %>%
  rowwise() %>%
  mutate(
    model_object = model_map[best_model],
    
    Adjusted_R2 = if (
      is.na(model_object) ||
      is.null(codon_results[[result]][[model_object]])
    ) {
      NA_real_
    } else {
      codon_results[[result]][[model_object]]$adj.r.squared
    },
    
    Position = codon,
    
    Codon_actual = if (codon == 1) {
      substr(pair, 1, 3)
    } else {
      substr(pair, 4, 6)
    }
  ) %>%
  ungroup() %>%
  select(
    Pair = pair,
    Position,
    Codon = Codon_actual,
    Best_model = best_model,
    Delta_AIC_to_GC3 = delta_to_gc,
    Adjusted_R2
  ) %>%
  mutate(
    Delta_AIC_to_GC3 = round(Delta_AIC_to_GC3, 2),
    Adjusted_R2 = round(Adjusted_R2, 3)
  ) %>%
  arrange(desc(Delta_AIC_to_GC3))

codon_appendix

write.csv(
  codon_appendix,
  "codon_appendix.csv",
  row.names = FALSE
)
