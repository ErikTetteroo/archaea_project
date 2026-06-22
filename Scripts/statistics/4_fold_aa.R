# four fold
# "Val" "Pro" "Thr" "Ala" "Gly"
v_aa <- aa_dat("Val")
table(v_aa$codon, v_aa$CM)

v_aa <- v_aa %>%
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "GUU" & CM == "M")
  ) %>%
  ungroup()

state <- v_aa %>%
  filter(codon == "GUG") %>%
  select(organism_id, state = CM)

dat <- v_aa %>%
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

ggplot(
  dat,
  aes(state, RSCU)
) +
  geom_boxplot() +
  facet_wrap(~ codon)



val_dat <- dat %>%
  mutate(
    GC_axis =
      (GUC + GUG) -
      (GUU + GUA),
    
    GUG_axis =
      GUG -
      (GUA + GUC + GUU)/3
  )

rownames(val_dat) <- val_dat$organism_id

fit_gc <- phylolm(
  GC_axis ~ GC3 + state,
  phy = tree,
  data = val_dat,
  model = "lambda"
)

fit_gug <- phylolm(
  GUG_axis ~ GC3 + state,
  phy = tree,
  data = val_dat,
  model = "lambda"
)


fit_gug1 <- phylolm(
  GUG_axis ~ GC3,
  phy = tree,
  data = val_dat,
  model = "lambda"
)


summary(fit_gc)
summary(fit_gug)

fit_gug$aic
fit_gug1$aic

fit_gc

ggplot(
  val_dat,
  aes(GC3, GUG_axis, color = state)
) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE)


# "Thr"
v_aa <- aa_dat("Thr")
table(v_aa$codon, v_aa$CM)

v_aa <- v_aa %>%
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
      !any(codon == "GUU" & CM == "M")
  ) %>%
  ungroup()

state <- v_aa %>%
  filter(codon == "ACG") %>%
  select(organism_id, state = CM)

dat <- v_aa %>%
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

ggplot(
  dat,
  aes(state, RSCU)
) +
  geom_boxplot() +
  facet_wrap(~ codon)



thr_dat <- dat %>%
  mutate(
    GC_axis =
      (ACC + ACG) -
      (ACU + ACA),
    
    ACG_axis =
      ACG -
      (ACC + ACU + ACA)/3
  )

rownames(thr_dat) <- thr_dat$organism_id

fit_gc <- phylolm(
  GC_axis ~ GC3 + state,
  phy = tree,
  data = thr_dat,
  model = "lambda"
)

fit_gug <- phylolm(
  ACG_axis ~ GC3 + state,
  phy = tree,
  data = thr_dat,
  model = "lambda"
)


fit_gug1 <- phylolm(
  ACG_axis ~ GC3,
  phy = tree,
  data = thr_dat,
  model = "lambda"
)


summary(fit_gc)
summary(fit_gug)

fit_gug$aic
fit_gug1$aic

fit_gc

ggplot(
  thr_dat,
  aes(GC3, ACG_axis, color = state)
) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE)

#"Ala"
v_aa <- aa_dat("Ala")
table(v_aa$codon, v_aa$CM)

v_aa <- v_aa %>%
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU")# &
      #!any(codon == "GUU" & CM == "M")
  ) %>%
  ungroup()

state <- v_aa %>%
  filter(codon == "GCG") %>%
  select(organism_id, state = CM)

dat <- v_aa %>%
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

ggplot(
  dat,
  aes(state, RSCU)
) +
  geom_boxplot() +
  facet_wrap(~ codon)



ala_dat <- dat %>%
  mutate(
    GC_axis =
      (GCC + GCG) -
      (GCU + GCA),
    
    GCG_axis =
      GCG -
      (GCC + GCU + GCA)/3
  )

rownames(ala_dat) <- ala_dat$organism_id

fit_gc <- phylolm(
  GC_axis ~ GC3 + state,
  phy = tree,
  data = ala_dat,
  model = "lambda"
)

fit_gug <- phylolm(
  GCG_axis ~ GC3 + state,
  phy = tree,
  data = ala_dat,
  model = "lambda"
)


fit_gug1 <- phylolm(
  GCG_axis ~ GC3,
  phy = tree,
  data = ala_dat,
  model = "lambda"
)


summary(fit_gc)
summary(fit_gug)

fit_gug$aic
fit_gug1$aic


ggplot(
  ala_dat,
  aes(GC3, GCG_axis, color = state)
) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE)

# "Gly"
v_aa <- aa_dat("Gly")
table(v_aa$codon, v_aa$CM)

v_aa <- v_aa %>%
  group_by(organism_id) %>%
  filter(
    !any(CM == "SU") &
    !any(codon == "GGU" & CM == "M")
  ) %>%
  ungroup()

state <- v_aa %>%
  filter(codon == "GGG") %>%
  select(organism_id, state = CM)

dat <- v_aa %>%
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

ggplot(
  dat,
  aes(state, RSCU)
) +
  geom_boxplot() +
  facet_wrap(~ codon)



gly_dat <- dat %>%
  mutate(
    GC_axis =
      (GGC + GGG) -
      (GGU + GGA),
    
    GGG_axis =
      GGG -
      (GGC + GGU + GGA)/3
  )

rownames(gly_dat) <- gly_dat$organism_id

fit_gc <- phylolm(
  GC_axis ~ GC3 + state,
  phy = tree,
  data = gly_dat,
  model = "lambda"
)

fit_gug <- phylolm(
  GGG_axis ~ GC3 + state,
  phy = tree,
  data = gly_dat,
  model = "lambda"
)


fit_gug1 <- phylolm(
  GGG_axis ~ GC3,
  phy = tree,
  data = gly_dat,
  model = "lambda"
)


summary(fit_gc)
summary(fit_gug)

fit_gug$aic
fit_gug1$aic


ggplot(
  gly_dat,
  aes(GC3, GGG_axis, color = state)
) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE)


results <- as.data.frame(coef(summary(fit_gug))["stateM", ])

ggplot(results,
       aes(
         y = "amino_acid",
         x = results[1,],
         xmin = results[1,] - 1.96*results[2,],
         xmax = results[1,] + 1.96*results[2,]
       )) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  geom_pointrange()

#"Pro"
v_aa <- aa_dat("Pro")
table(v_aa$codon, v_aa$CM)

v_aa <- v_aa %>%
  group_by(organism_id) %>%
  filter(
    #!any(CM == "SU")# &
    !any(codon == "CCU" & CM == "M")
  ) %>%
  ungroup()

state <- v_aa %>%
  filter(codon == "GCG") %>%
  select(organism_id, state = CM)

pro_state1 <- v_aa %>%
  filter(codon == "CCG") %>%
  select(organism_id,
         CCG_state = CM)

pro_state2 <- v_aa %>%
  filter(codon == "CCC") %>%
  select(organism_id,
         CCC_state = CM)

pro_state3 <- v_aa %>%
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

pro_dat <- v_aa %>%
  left_join(pro_states, by = "organism_id")

table(pro_dat$CCG_state, pro_dat$SU_state)

ggplot(pro_states, aes(GC3)) +
  geom_density(aes(fill = CCG_state), alpha = 0.4)

dat <- v_aa %>%
  inner_join(pro_states, by = "organism_id") %>%
  select(
    organism_id,
    codon,
    RSCU,
    GC3,
    CCG_state,
    CCC_state,
    CCU_state
  ) %>%
  pivot_wider(
    names_from = codon,
    values_from = RSCU
  )

ggplot(
  dat,
  aes(CCU_state, RSCU)
) +
  geom_boxplot() +
  facet_wrap(~ codon)

ggplot(dat, aes(GC3)) +
  geom_density(aes(fill = CCC_state), alpha = 0.4)

dat <- dat %>%
  mutate(
    pro_state = case_when(
      CCG_state == "M" & CCC_state == "M" ~ "high",
      CCG_state == "GU" & CCC_state == "M" ~ "mid",
      CCG_state == "GU" & CCC_state == "SU" ~ "low",
      TRUE ~ NA_character_
    )
  )

dat$pro_state <- factor(
  dat$pro_state,
  levels = c("low", "mid", "high"),
  ordered = TRUE
)

pro_dat <- dat %>%
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

fit_gc <- phylolm(
  GC_axis ~ GC3 + CCG_state,
  phy = tree,
  data = pro_dat,
  model = "lambda"
)

fit_gug <- phylolm(
  CCG_axis ~ GC3 + CCG_state,
  phy = tree,
  data = pro_dat,
  model = "lambda"
)


fit_gug1 <- phylolm(
  CCG_axis ~ GC3,
  phy = tree,
  data = pro_dat,
  model = "lambda"
)

fit_su <- phylolm(
  SU_axis ~ GC3 + pro_state,
  phy = tree,
  data = pro_dat,
  model = "lambda"
)

fit_su1 <- phylolm(
  SU_axis ~ GC3,
  phy = tree,
  data = pro_dat,
  model = "lambda"
)


summary(fit_gc)
summary(fit_gug)

fit_gug$aic
fit_gug1$aic

summary(fit_su)

fit_su$aic
fit_su1$aic

ggplot(
  pro_dat,
  aes(GC3, SU_axis, color = pro_state)
) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE)

ggplot(pro_dat, aes(GC3, fill = pro_state)) +
  geom_density(alpha = 0.4)
