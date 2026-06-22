# subsets the codon usage data of codons from a certain amino acid
# filters out organisms with codons that aren't covered
aa_subset <- function(data,aa) {
  merged %>%
    filter(
      amino_acid == aa
    ) %>%
    group_by(organism_id) %>%
    filter(!any(is.na(CM))) %>%
    ungroup()
}

# two-fold analysis

analyze_twofold <- function(aa_ss, variable_codon, tree) {
  
  state <- aa_ss %>%
    filter(codon == variable_codon) %>%
    select(organism_id, state = CM)
  
  dat <- aa_ss %>%
    left_join(state, by = "organism_id") %>%
    select(
      organism_id,
      codon,
      RSCU,
      GC3,
      state
    ) %>%
    pivot_wider(
      names_from = codon,
      values_from = RSCU
    )
  
  codons <- setdiff(names(dat),
                    c("organism_id",
                      "GC3",
                      "state"))
  
  ifelse(codons[1]==variable_codon,
         dat$delta <- dat[[codons[1]]] -
           dat[[codons[2]]],
         dat$delta <- dat[[codons[2]]] -
           dat[[codons[1]]])
  
  rownames(dat) <- dat$organism_id
  
  fit_gc <- phylolm(
    delta ~ GC3,
    data = dat,
    phy = tree,
    model = "lambda"
  )
  
  fit_quad <- phylolm(
    delta ~ GC3 + I(GC3^2),
    data = dat,
    phy = tree,
    model = "lambda"
  )
  
  fit_full <- phylolm(
    delta ~ GC3 + state,
    data = dat,
    phy = tree,
    model = "lambda"
  )
  
  fit_quad_full <- phylolm(
    delta ~ GC3 + I(GC3^2) + state,
    data = dat,
    phy = tree,
    model = "lambda"
  )
  
  list(
    dat = dat,
    fit_gc = fit_gc,
    fit_quad = fit_quad,
    fit_full = fit_full,
    fit_quad_full = fit_quad_full
  )
}

best_model <- function(res) {
  
  aics <- c(
    gc = res$fit_gc$aic,
    quad = res$fit_quad$aic,
    full = res$fit_full$aic,
    quad_full = res$fit_quad_full$aic
  )
  
  names(which.min(aics))
}

model_table <- function(res) {
  
  data.frame(
    model = c(
      "GC3",
      "GC3 + GC3²",
      "GC3 + state",
      "GC3 + GC3² + state"
    ),
    AIC = c(
      res$fit_gc$aic,
      res$fit_quad$aic,
      res$fit_full$aic,
      res$fit_quad_full$aic
    )
  )
}

coef_table <- function(model) {
  
  as.data.frame(summary(model)$coefficients) %>%
    tibble::rownames_to_column("term")
}

save_twofold <- function(res, aa) {
  
  dir.create(aa, showWarnings = FALSE)
  
  #---------------------------
  # choose best model
  #---------------------------
  
  best <- switch(
    best_model(res),
    gc = res$fit_gc,
    quad = res$fit_quad,
    full = res$fit_full,
    quad_full = res$fit_quad_full
  )
  
  #---------------------------
  # model selection
  #---------------------------
  
  write_csv(
    model_table(res),
    file.path(aa, "model_selection.csv")
  )
  
  #---------------------------
  # model summary
  #---------------------------
  
  capture.output(
    summary(best),
    file = file.path(aa, "model_summary.txt")
  )
  
  #---------------------------
  # predictions
  #---------------------------
  
  pred_dat <- expand.grid(
    GC3 = seq(
      min(res$dat$GC3),
      max(res$dat$GC3),
      length.out = 200
    ),
    state = c("GU", "M")
  )
  
  b <- coef(res$fit_quad_full)
  
  pred_dat <- pred_dat %>%
    mutate(
      pred =
        b["(Intercept)"] +
        b["GC3"] * GC3 +
        b["I(GC3^2)"] * GC3^2 +
        ifelse(state == "M",
               b["stateM"],
               0)
    )
  
  V <- vcov(res$fit_quad_full)
  
  X <- model.matrix(
    ~ GC3 + I(GC3^2) + state,
    pred_dat
  )
  
  pred_dat$se <- sqrt(
    diag(X %*% V %*% t(X))
  )
  
  pred_dat <- pred_dat %>%
    mutate(
      lower = pred - 1.96 * se,
      upper = pred + 1.96 * se
    )
  
  #---------------------------
  # GC3 plot
  #---------------------------
  
  p1 <- ggplot(res$dat,
               aes(GC3,
                   delta,
                   color = state)) +
    
    geom_point(alpha = 0.3) +
    
    geom_ribbon(
      data = pred_dat,
      aes(
        x = GC3,
        ymin = lower,
        ymax = upper,
        fill = state
      ),
      alpha = 0.2,
      colour = NA,
      inherit.aes = FALSE
    ) +
    
    geom_line(
      data = pred_dat,
      aes(
        x = GC3,
        y = pred,
        color = state
      ),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    
    coord_cartesian(
      ylim = c(-2, 2)
    )
  
  ggsave(
    file.path(aa, "GC3_vs_delta.pdf"),
    p1,
    width = 6,
    height = 4
  )
  
  #---------------------------
  # residuals
  #---------------------------
  
  pdf(file.path(aa, "residuals.pdf"))
  
  plot(
    fitted(best),
    residuals(best),
    xlab = "Fitted",
    ylab = "Residuals"
  )
  
  abline(h = 0, lty = 2)
  
  dev.off()
  
  #---------------------------
  # QQ plot
  #---------------------------
  
  pdf(file.path(aa, "qqplot.pdf"))
  
  qqnorm(residuals(best))
  qqline(residuals(best))
  
  dev.off()
  
}

summary_table <- function(res, aa) {
  
  best <- switch(
    best_model(res),
    gc = res$fit_gc,
    quad = res$fit_quad,
    full = res$fit_full
  )
  
  data.frame(
    AA = aa,
    lambda = best$optpar,
    AIC_gc = res$fit_gc$aic,
    AIC_quad = res$fit_quad$aic,
    AIC_full = res$fit_full$aic,
    AIC_quad_full = res$fit_quad_full$aic
  )
}

# four_fold analysis

fourfold_prep <- function(aa_ss, variable_codon) {
  states <- aa_ss %>%
    filter(codon == variable_codon) %>%
    select(organism_id, state = CM)
  
  aa_ss %>%
    left_join(states, by = "organism_id") %>%
    select(
      organism_id,
      codon,
      RSCU,
      GC3,
      state
    ) 
}

state_rscu_boxplot <- function(aa_ss, title) {
ggplot(
  aa_ss,
  aes(state, RSCU)
) +
  geom_boxplot() +
  facet_wrap(~ codon) +
  ggtitle(title) +
  coord_cartesian(ylim = c(0, 4)) 
}

four_fold_axes <- function(aa_ss, variable_codon) {
  
  # remove the third-position base
  prefix <- substr(variable_codon, 1, 2)
  
  codons <- c(
    U = paste0(prefix, "U"),
    C = paste0(prefix, "C"),
    A = paste0(prefix, "A"),
    G = paste0(prefix, "G")
  )
  
  others <- setdiff(codons, variable_codon)
  
  aa_ss %>%
    mutate(
      GC_axis =
        .data[[codons["C"]]] +
        .data[[codons["G"]]] -
        .data[[codons["U"]]] -
        .data[[codons["A"]]],
      
      var_axis =
        .data[[variable_codon]] -
        rowMeans(across(all_of(others)))
    )
}

analyze_fourfold <- function(dat, tree) {
  
  fit_gc1 <- phylolm(
    GC_axis ~ GC3,
    data = dat,
    phy = tree,
    model = "lambda"
  )
  
  fit_gc2 <- phylolm(
    GC_axis ~ GC3 + I(GC3^2),
    data = dat,
    phy = tree,
    model = "lambda"
  )
  
  fit_gc3 <- phylolm(
    GC_axis ~ GC3 + I(GC3^2) + state,
    data = dat,
    phy = tree,
    model = "lambda"
  )
  
  fit_var1 <- phylolm(var_axis ~ GC3,
                      data = dat,
                      phy = tree,
                      model = "lambda")
  
  fit_var2 <- phylolm(var_axis ~ GC3 + I(GC3^2),
                      data = dat,
                      phy = tree,
                      model = "lambda")
  
  fit_var3 <- phylolm(var_axis ~ GC3 + state,
                      data = dat,
                      phy = tree,
                      model = "lambda")
  
  fit_var4 <- phylolm(var_axis ~ GC3 + I(GC3^2) + state,
                      data = dat,
                      phy = tree,
                      model = "lambda")
  
  list(
    dat = dat,
    
    fit_gc1 = fit_gc1,
    fit_gc2 = fit_gc2,
    fit_gc3 = fit_gc3,
    
    fit_var1 = fit_var1,
    fit_var2 = fit_var2,
    fit_var3 = fit_var3,
    fit_var4 = fit_var4
  )
}

delta_aic <- function(x) {
  x - min(x)
}

best_fourfold <- function(res) {
  
  gc_aic <- c(
    linear = res$fit_gc1$aic,
    quadratic = res$fit_gc2$aic,
    quadratic_state = res$fit_gc3$aic
  )
  
  var_aic <- c(
    linear = res$fit_var1$aic,
    quadratic = res$fit_var2$aic,
    linear_state = res$fit_var3$aic,
    quadratic_state = res$fit_var4$aic
  )
  
  best_gc <- names(which.min(gc_aic))
  best_var <- names(which.min(var_aic))
  
  list(
    gc_model = best_gc,
    gc_fit = res[[c(
      linear = "fit_gc1",
      quadratic = "fit_gc2",
      quadratic_state = "fit_gc3"
    )[best_gc]]],
    
    var_model = best_var,
    var_fit = res[[c(
      linear = "fit_var1",
      quadratic = "fit_var2",
      linear_state = "fit_var3",
      quadratic_state = "fit_var4"
    )[best_var]]],
    
    gc_aic = gc_aic,
    var_aic = var_aic
  )
}

save_fourfold <- function(fit,
                          best,
                          aa) {
  
  dir.create(aa, showWarnings = FALSE)
  
  ## ----------------------------
  ## AIC table
  ## ----------------------------
  
  aic_table <- bind_rows(
    tibble(
      axis = "GC_axis",
      model = names(best$gc_aic),
      AIC = unname(best$gc_aic),
      deltaAIC = unname(best$gc_delta)
    ),
    tibble(
      axis = "variable_axis",
      model = names(best$var_aic),
      AIC = unname(best$var_aic),
      deltaAIC = unname(best$var_delta)
    )
  )
  
  write_csv(
    aic_table,
    file.path(aa, "model_selection.csv")
  )
  
  ## ----------------------------
  ## Model summaries
  ## ----------------------------
  
  sink(file.path(aa, "model_summary.txt"))
  
  cat("GC AXIS MODEL\n")
  cat("====================\n\n")
  print(summary(best$gc_fit))
  
  cat("\n\n")
  cat("VARIABLE AXIS MODEL\n")
  cat("====================\n\n")
  print(summary(best$var_fit))
  
  sink()
  
  ## ----------------------------
  ## Variable axis plot
  ## ----------------------------
  
  pred_dat <- expand.grid(
    GC3 = seq(
      min(fit$dat$GC3),
      max(fit$dat$GC3),
      length.out = 200
    ),
    state = unique(fit$dat$state)
  )
  
  fit_vis <- fit$fit_var4
  
  X <- model.matrix(
    ~ GC3 + I(GC3^2) + state,
    pred_dat
  )
  
  b <- coef(fit_vis)
  
  pred_dat$pred <- X %*% b
  
  V <- vcov(fit_vis)
  
  pred_dat$se <- sqrt(
    diag(X %*% V %*% t(X))
  )
  
  pred_dat <- pred_dat %>%
    mutate(
      lower = pred - 1.96 * se,
      upper = pred + 1.96 * se
    )
  
  pdf(file.path(aa, "variable_axis.pdf"),
      width = 6,
      height = 5)
  
  p <- ggplot(
    fit$dat,
    aes(
      GC3,
      var_axis,
      color = state
    )
  ) +
    geom_point(alpha = 0.3) +
    
    geom_ribbon(
      data = pred_dat,
      aes(
        x = GC3,
        ymin = lower,
        ymax = upper,
        fill = state
      ),
      alpha = 0.2,
      colour = NA,
      inherit.aes = FALSE
    ) +
    
    geom_line(
      data = pred_dat,
      aes(
        x = GC3,
        y = pred,
        color = state
      ),
      linewidth = 1,
      inherit.aes = FALSE
    )
  print(p)
  
  dev.off()
  
  ## ----------------------------
  ## Residual plot
  ## ----------------------------
  
  pdf(file.path(aa, "residuals.pdf"),
      width = 5,
      height = 5)
  
  plot(
    fitted(best$var_fit),
    residuals(best$var_fit),
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  
  abline(
    h = 0,
    lty = 2
  )
  
  dev.off()
  
  ## ----------------------------
  ## QQ plot
  ## ----------------------------
  
  pdf(file.path(aa, "qqplot.pdf"),
      width = 5,
      height = 5)
  
  qqnorm(residuals(best$var_fit))
  qqline(residuals(best$var_fit))
  
  dev.off()
}
