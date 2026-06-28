# load project
source("R/load_project.R")

# read in files
merged <- read_csv("Data/cleaned_data/codon_usage_data_stats_ready.csv")
tree <- read.tree("Data/cleaned_data/cleaned_archaea.tree")
valid_amino_acids <- read_csv("Data/cleaned_data/aa_to_analyse.csv")

# subset the interesting 4 fold amino acids
four_fold <- valid_amino_acids$aa[valid_amino_acids$Freq==4]

four_fold_subsets <- list()

for (n in four_fold) {
  four_fold_subsets[[n]] <- aa_subset(merged, n)
}

#-----------------------------------------------------------------
# Check which codons have sufficient variation and filter outliers
#-----------------------------------------------------------------
# Val
table(four_fold_subsets$Val$codon, four_fold_subsets$Val$CM) # GUG varies with 36 alternative samples
four_fold_subsets$Val <- four_fold_subsets$Val %>%           # filter super wobble cases and GUU matches
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "GUU" & CM == "M")
  ) %>%
  ungroup()
#Pro
table(four_fold_subsets$Pro$codon, four_fold_subsets$Pro$CM) # CCG varies with 74 alternative samples
four_fold_subsets$Pro <- four_fold_subsets$Pro %>%           # Super wobble variation is present in 25 samples
  group_by(organism_id) %>%                                  # filter the single CCU match  
  filter(
    !any(codon == "CCU" & CM == "M")
  ) %>%
  ungroup()
#Thr
table(four_fold_subsets$Thr$codon, four_fold_subsets$Thr$CM) # ACG varies with 29 alternative samples
#Ala
table(four_fold_subsets$Ala$codon, four_fold_subsets$Ala$CM) # GCG varies with 47 alternative samples
four_fold_subsets$Ala <- four_fold_subsets$Ala %>%           # filter super wobble cases  
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU")
  ) %>%
  ungroup()
#Gly
table(four_fold_subsets$Gly$codon, four_fold_subsets$Gly$CM) # GGG varies with 84 alternative samples
four_fold_subsets$Gly <- four_fold_subsets$Gly %>%           #filter super wobble cases and GGU matches
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "GGU" & CM == "M")
  ) %>%
  ungroup()

#-----------------------------------------------------------------
# split organisms based on the state of their varying codons
#-----------------------------------------------------------------
four_fold_states <- list()

four_fold_states[["Val"]] <- fourfold_prep(four_fold_subsets$Val,"GUG")
four_fold_states[["Thr"]] <- fourfold_prep(four_fold_subsets$Thr,"ACG")
four_fold_states[["Ala"]] <- fourfold_prep(four_fold_subsets$Ala,"GCG")
four_fold_states[["Gly"]] <- fourfold_prep(four_fold_subsets$Gly,"GGG")

state_rscu_boxplot(four_fold_states$Val, "Val")
state_rscu_boxplot(four_fold_states$Thr, "Thr")
state_rscu_boxplot(four_fold_states$Ala, "Ala")
state_rscu_boxplot(four_fold_states$Gly, "Gly")

for (i in 1:length(four_fold_states)) {
  four_fold_states[[i]] <- four_fold_states[[i]] %>%
    pivot_wider(
      names_from = codon,
      values_from = RSCU
    )
}

#-----------------------------------------------------------------
# split respons variable in 2 axes, CG vs AU; and variable codon vs the rest
#-----------------------------------------------------------------
four_fold_ready <- list()

four_fold_ready[["Val"]] <- four_fold_axes(four_fold_states$Val,"GUG")
four_fold_ready[["Thr"]] <- four_fold_axes(four_fold_states$Thr,"ACG")
four_fold_ready[["Ala"]] <- four_fold_axes(four_fold_states$Ala,"GCG")
four_fold_ready[["Gly"]] <- four_fold_axes(four_fold_states$Gly,"GGG")

for (i in 1:length(four_fold_ready)) {
  rownames(four_fold_ready[[i]]) <- four_fold_ready[[i]]$organism_id
}

#-----------------------------------------------------------------
# fit the models
#-----------------------------------------------------------------

# Val
fit_val <- analyze_fourfold(four_fold_ready$Val,tree)
best_val <- best_fourfold(fit_val)
best_val$gc_delta <-
  delta_aic(best_val$gc_aic)

best_val$var_delta <-
  delta_aic(best_val$var_aic)

# Thr
fit_thr <- analyze_fourfold(four_fold_ready$Thr,tree)
best_thr <- best_fourfold(fit_thr)
best_thr$gc_delta <-
  delta_aic(best_thr$gc_aic)

best_thr$var_delta <-
  delta_aic(best_thr$var_aic)

# Ala
fit_ala <- analyze_fourfold(four_fold_ready$Ala,tree)
best_ala <- best_fourfold(fit_ala)
best_ala$gc_delta <-
  delta_aic(best_ala$gc_aic)

best_ala$var_delta <-
  delta_aic(best_ala$var_aic)

# Gly
fit_gly <- analyze_fourfold(four_fold_ready$Gly,tree)
best_gly <- best_fourfold(fit_gly)
best_gly$gc_delta <-
  delta_aic(best_gly$gc_aic)

best_gly$var_delta <-
  delta_aic(best_gly$var_aic)

# save outputs
save_fourfold(
  fit = fit_val,
  best = best_val,
  aa = "Plots/aa/four_fold/Val"
)
save_fourfold(
  fit = fit_thr,
  best = best_thr,
  aa = "Plots/aa/four_fold/Thr"
)
save_fourfold(
  fit = fit_ala,
  best = best_ala,
  aa = "Plots/aa/four_fold/Ala"
)
save_fourfold(
  fit = fit_gly,
  best = best_gly,
  aa = "Plots/aa/four_fold/Gly"
)

#---------------------------------------------------------------
# Proline
#---------------------------------------------------------------
pro_state1 <- four_fold_subsets$Pro %>%
  filter(codon == "CCG") %>%
  select(organism_id,
         CCG_state = CM)

pro_state2 <- four_fold_subsets$Pro %>%
  filter(codon == "CCC") %>%
  select(organism_id,
         CCC_state = CM)

pro_state3 <- four_fold_subsets$Pro %>%
  filter(codon == "CCU") %>%
  select(organism_id,
         CCU_state = CM)

pro_states <- pro_state1 %>%
  left_join(
    pro_state2,
    by = "organism_id"
  ) %>%
  left_join(
    pro_state3,
    by = "organism_id"
  )

pro_states %>%
  count(CCG_state, CCC_state, CCU_state)

# filter the 2 organisms that have M CCG and SU
pro_states <- pro_states %>%
  group_by(organism_id) %>%
  filter(
    !any(CCG_state == "M" & CCC_state == "SU")
  ) %>%
  ungroup()

pro_states <- pro_states %>%
  mutate(
    CCG_state = factor(CCG_state),
    SU_state  = ifelse(CCC_state == "SU", "SU", "M"),
    SU_state  = factor(SU_state)
  )

pro_dat <- four_fold_subsets$Pro %>%
  inner_join(pro_states, by = "organism_id") %>%
  select(
    organism_id,
    codon,
    RSCU,
    GC3,
    CCG_state,
    SU_state
  ) 

table(pro_dat$CCG_state, pro_dat$SU_state)


pro_dat <- pro_dat %>%
  mutate(
    pro_state = case_when(
      CCG_state == "M" ~ "M",
      CCG_state == "GU" & SU_state == "M" ~ "GU",
      CCG_state == "GU" & SU_state == "SU" ~ "GU_SU"
    ),
    pro_state = factor(pro_state)
  ) 

pro_dat <- pro_dat %>%
  pivot_wider(
    names_from = codon,
    values_from = RSCU
  )

pro_dat <- pro_dat %>%
  mutate(
    GC_axis =
      (CCC + CCG) -
      (CCU + CCA),
    
    CCG_axis =
      CCG -
      (CCC + CCU + CCA)/3,
    
    SU_axis =
      (CCC + CCU) -
      (CCG + CCA)
  )

rownames(pro_dat) <- pro_dat$organism_id

fit_gc1 <- phylolm(
  GC_axis ~ GC3,
  data = pro_dat,
  phy = tree,
  model = "lambda"
)

fit_gc2 <- phylolm(
  GC_axis ~ GC3 + I(GC3^2),
  data = pro_dat,
  phy = tree,
  model = "lambda"
)

fit_gc3 <- phylolm(
  GC_axis ~ GC3 + I(GC3^2) + pro_state,
  data = pro_dat,
  phy = tree,
  model = "lambda"
)

fit_var1 <- phylolm(CCG_axis ~ GC3,
                    data = pro_dat,
                    phy = tree,
                    model = "lambda")

fit_var2 <- phylolm(CCG_axis ~ GC3 + I(GC3^2),
                    data = pro_dat,
                    phy = tree,
                    model = "lambda")

fit_var3 <- phylolm(CCG_axis ~ GC3 + pro_state,
                    data = pro_dat,
                    phy = tree,
                    model = "lambda")

fit_var4 <- phylolm(CCG_axis ~ GC3 + I(GC3^2) + pro_state,
                    data = pro_dat,
                    phy = tree,
                    model = "lambda")

pro_fit <- list(
  dat = pro_dat,
  
  fit_gc1 = fit_gc1,
  fit_gc2 = fit_gc2,
  fit_gc3 = fit_gc3,
  
  fit_var1 = fit_var1,
  fit_var2 = fit_var2,
  fit_var3 = fit_var3,
  fit_var4 = fit_var4
)

best_pro <- best_fourfold(pro_fit)
best_pro$gc_delta <-
  delta_aic(best_pro$gc_aic)

best_pro$var_delta <-
  delta_aic(best_pro$var_aic)

## ----------------------------
## Variable axis plot
## ----------------------------

pred_dat <- expand.grid(
  GC3 = seq(
    min(pro_fit$dat$GC3),
    max(pro_fit$dat$GC3),
    length.out = 200
  ),
  state = unique(pro_fit$dat$pro_state)
)

fit_vis <- pro_fit$fit_var3

X <- model.matrix(
  ~ GC3 + state,
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

jpeg(file.path("Plots/aa/four_fold/Pro/variable_axis.jpeg"),
    width = 600,
    height = 500)

p <- ggplot(
  pro_fit$dat,
  aes(
    GC3,
    CCG_axis,
    color = pro_state
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

pdf(file.path("Plots/aa/four_fold/Pro/residuals.pdf"),
    width = 5,
    height = 5)

plot(
  fitted(best_pro$var_fit),
  residuals(best_pro$var_fit),
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

pdf(file.path("Plots/aa/four_fold/Pro/qqplot.pdf"),
    width = 5,
    height = 5)

qqnorm(residuals(best_pro$var_fit))
qqline(residuals(best_pro$var_fit))

dev.off()


