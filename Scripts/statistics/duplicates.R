# load project
source("R/load_project.R")

# read in files
merged <- read_csv("Data/cleaned_data/codon_usage_data_stats_ready.csv")
tree <- read.tree("Data/cleaned_data/cleaned_archaea.tree")

# exploration
dup_summary <-
  merged %>%
  group_by(codon) %>%
  summarise(
    n_states = n_distinct(M),
    min_M = min(M),
    max_M = max(M),
    mean_M = mean(M),
    sd_M = sd(M)
  ) %>%
  arrange(desc(max_M))

dup_summary <- dup_summary[1:10,]

dup_subset <- merged[merged$codon %in% dup_var_codons,]
dup_subset_a <- merged[merged$amino_acid %in% dup_subset$amino_acid,]



ggplot(
  dup_summary,
  aes(
    reorder(codon, max_M),
    max_M
  )
) +
  geom_col() +
  coord_flip()

ggplot(
  dup_subset,
  aes(factor(M))
) +
  geom_bar() +
  facet_wrap(~codon)

test <- dup_subset[dup_subset$codon=="UUA",]

for (n in unique(dup_subset$codon)) {
  print(summary(dup_subset$M[dup_subset$codon==n]))
  
}

t <- dup_subset %>%
  count(codon, M) 

#"GCA"

dup_var_codons <- c("UGC","GAC")



make_dup_twofold <- function(dat, variable_codon) {
  
  other_codon <-
    setdiff(unique(dat$codon), variable_codon)
  
  dup_info <-
    dat %>%
    filter(codon == variable_codon) %>%
    select(
      organism_id,
      M
    )
  
  dat %>%
    select(
      organism_id,
      codon,
      RSCU,
      GC3
    ) %>%
    pivot_wider(
      names_from = codon,
      values_from = RSCU
    ) %>%
    left_join(
      dup_info,
      by = "organism_id"
    ) %>%
    mutate(
      delta =
        .data[[variable_codon]] -
        .data[[other_codon]]
    ) %>%
    column_to_rownames("organism_id")
}

asp_subset <-
  dup_subset_a %>%
  filter(
    amino_acid == "Asp"
  ) %>%
  group_by(organism_id) %>%
  filter(
    M[codon == "GAC"] > 0
  ) %>%
  ungroup()

cys_subset <-
  dup_subset_a %>%
  filter(
    amino_acid == "Cys"
  ) %>%
  group_by(organism_id) %>%
  filter(
    M[codon == "UGC"] > 0
  ) %>%
  ungroup()

asp_dat <- make_dup_twofold(
  asp_subset,
  "GAC"
)

cys_dat <- make_dup_twofold(
  cys_subset,
  "UGC"
)


ggplot(
  asp_dat,
  aes(factor(M), GC3)
) +
  geom_boxplot()


ggplot(
  cys_dat,
  aes(factor(M), GC3)
) +
  geom_boxplot()

ggplot(
  cys_dat,
  aes(factor(M), delta)
) +
  geom_boxplot()

ggplot(
  cys_dat,
  aes(GC3, delta, color = factor(M))
) +
  geom_point(alpha = 0.5)

cor.test(asp_dat$GC3, asp_dat$M, method = "spearman")

cor.test(cys_dat$GC3, cys_dat$M, method = "spearman")

ggplot(
  cys_dat,
  aes(M, delta)
) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  geom_boxplot(
    aes(group = M),
    alpha = 0.3
  )

ggplot(
  cys_dat,
  aes(GC3, delta, color = factor(M))
) +
  geom_point(alpha = 0.5)


fit1 <- phylolm(
  delta ~ GC3,
  data = asp_dat,
  phy = tree,
  model = "lambda"
)

fit2 <- phylolm(
  delta ~ GC3 + I(GC3^2),
  data = asp_dat,
  phy = tree,
  model = "lambda"
)

fit3 <- phylolm(
  delta ~ GC3 + M,
  data = asp_dat,
  phy = tree,
  model = "lambda"
)

fit3l <- phylolm(
  delta ~ GC3 + log(M),
  data = asp_dat,
  phy = tree,
  model = "lambda"
)

fit4 <- phylolm(
  delta ~ GC3 + I(GC3^2) + M,
  data = asp_dat,
  phy = tree,
  model = "lambda"
)

fit4l <- phylolm(
  delta ~ GC3 + I(GC3^2) + log(M),
  data = asp_dat,
  phy = tree,
  model = "lambda"
)

fit1$aic
fit2$aic
fit3$aic
fit3l$aic
fit4$aic
fit4l$aic

summary(fit3l)

fit1 <- phylolm(
  delta ~ GC3,
  data = cys_dat,
  phy = tree,
  model = "lambda"
)

fit2 <- phylolm(
  delta ~ GC3 + I(GC3^2),
  data = cys_dat,
  phy = tree,
  model = "lambda"
)

fit3 <- phylolm(
  delta ~ GC3 + M,
  data = cys_dat,
  phy = tree,
  model = "lambda"
)

fit3l <- phylolm(
  delta ~ GC3 + log(M),
  data = cys_dat,
  phy = tree,
  model = "lambda"
)

fit4 <- phylolm(
  delta ~ GC3 + I(GC3^2) + M,
  data = cys_dat,
  phy = tree,
  model = "lambda"
)

fit4l <- phylolm(
  delta ~ GC3 + I(GC3^2) + log(M),
  data = cys_dat,
  phy = tree,
  model = "lambda"
)

fit1$aic
fit2$aic
fit3$aic
fit3l$aic
fit4$aic
fit4l$aic

summary(fit4l)


sd(log(cys_dat$M))
sd(cys_dat$GC3)

beta_M  <- coef(fit4l)["log(M)"] * sd(log(cys_dat$M))
beta_GC <- coef(fit4l)["GC3"] * sd(cys_dat$GC3)

#-----------------------------------------------------------------

pred_cys <- expand.grid(
  GC3 = seq(
    min(cys_dat$GC3),
    max(cys_dat$GC3),
    length.out = 200
  ),
  M = c(1,2,4)
)

X <- model.matrix(
  ~ GC3 + I(GC3^2) + log(M),
  pred_cys
)

b <- coef(fit4l)

pred_cys$pred <- X %*% b

V <- vcov(fit4l)

pred_cys$se <- sqrt(
  diag(X %*% V %*% t(X))
)

pred_cys <- pred_cys %>%
  mutate(
    lower = pred - 1.96*se,
    upper = pred + 1.96*se,
    M = factor(M)
  )

ggplot(
  cys_dat,
  aes(GC3, delta)
) +
  
  geom_point(
    aes(color = factor(M)),
    alpha = 0.3
  ) +
  
  geom_ribbon(
    data = pred_cys,
    aes(
      x = GC3,
      ymin = lower,
      ymax = upper,
      fill = M,
      group = M
    ),
    alpha = 0.2,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = pred_cys,
    aes(
      x = GC3,
      y = pred,
      color = M,
      group = M
    ),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  
  coord_cartesian(
    ylim = c(-2, 2)
  ) +
  
  labs(
    x = "GC3",
    y = expression(delta),
    color = "Copies",
    fill = "Copies"
  ) +
  
  theme_classic()

#-------------------------------------------------------------

pred_asp <- expand.grid(
  GC3 = seq(
    min(asp_dat$GC3),
    max(asp_dat$GC3),
    length.out = 200
  ),
  M = c(1,2,4)
)

X <- model.matrix(
  ~ GC3 + log(M),
  pred_asp
)

b <- coef(fit3l)

pred_asp$pred <- X %*% b

V <- vcov(fit3l)

pred_asp$se <- sqrt(
  diag(X %*% V %*% t(X))
)

pred_asp <- pred_asp %>%
  mutate(
    lower = pred - 1.96*se,
    upper = pred + 1.96*se,
    M = factor(M)
  )

ggplot(
  asp_dat,
  aes(GC3, delta)
) +
  
  geom_point(
    aes(color = factor(M)),
    alpha = 0.3
  ) +
  
  geom_ribbon(
    data = pred_asp,
    aes(
      x = GC3,
      ymin = lower,
      ymax = upper,
      fill = M,
      group = M
    ),
    alpha = 0.2,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = pred_asp,
    aes(
      x = GC3,
      y = pred,
      color = M,
      group = M
    ),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  
  coord_cartesian(
    ylim = c(-2, 2)
  ) +
  
  labs(
    x = "GC3",
    y = expression(delta),
    color = "Copies",
    fill = "Copies"
  ) +
  
  theme_classic()
