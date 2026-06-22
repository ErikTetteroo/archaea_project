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
  
  ############################
  # choose best model
  ############################
  
  best <- switch(
    best_model(res),
    gc = res$fit_gc,
    quad = res$fit_quad,
    full = res$fit_full
  )
  
  ############################
  # model selection
  ############################
  
  write_csv(
    model_table(res),
    file.path(aa, "model_selection.csv")
  )
  
  ############################
  # model summary
  ############################
  
  capture.output(
    summary(best),
    file = file.path(aa, "model_summary.txt")
  )
  
  ############################
  # GC3 plot
  ############################
  
  p1 <- ggplot(
    res$dat,
    aes(GC3, delta, color = state)
  ) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess")
  
  ggsave(
    file.path(aa, "GC3_vs_delta.pdf"),
    p1,
    width = 6,
    height = 4
  )
  
  ############################
  # residuals
  ############################
  
  pdf(file.path(aa, "residuals.pdf"))
  
  plot(
    fitted(best),
    residuals(best),
    xlab = "Fitted",
    ylab = "Residuals"
  )
  
  abline(h = 0, lty = 2)
  
  dev.off()
  
  ############################
  # QQ plot
  ############################
  
  pdf(file.path(aa, "qqplot.pdf"))
  
  qqnorm(residuals(best))
  qqline(residuals(best))
  
  dev.off()
  
}