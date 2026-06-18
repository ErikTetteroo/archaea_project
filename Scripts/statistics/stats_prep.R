# load project
source("R/load_project.R")

merged <- read_csv("Data/cleaned_data/merged_codon_usage_data_with_pcoa_taxonomy.csv")

# identify interesting codons 
merged |>
  group_by(codon) |>
  summarise(
    n_M  = sum(CM == "M", na.rm = TRUE),
    n_GU = sum(CM == "GU", na.rm = TRUE),
    n_SU = sum(CM == "SU", na.rm = TRUE),
    n_M2 = sum(CM == "M2", na.rm = TRUE)
  ) |>
  print(n=63)

valid_codons <- merged %>%
  filter(!is.na(CM)) %>%
  count(codon, CM) %>%
  group_by(codon) %>%
  summarise(
    n_categories_20 = sum(n >= 20)
  ) %>%
  filter(n_categories_20 >= 2) %>%
  pull(codon)

table <- give_codon_table()
valid_amino_acids <- unique(table$aa[table$codon_rna %in% valid_codons])

#plot
merged_c <- merged[!is.na(merged$CM),]

cm_prop <- merged_c |>
  summarise(
    n = n(),
    .by = c(codon, CM)
  ) |>
  group_by(codon) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

ggplot(cm_prop, aes(codon, prop, fill = CM)) +
  geom_col() +
  labs(
    x = "Codon",
    y = "Proportion",
    fill = "Coverage type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )


aa_dat <- function(aa) {
  merged |>
    filter(
      in_tree,
      amino_acid == aa
    ) |>
    group_by(organism_id) |>
    filter(!any(is.na(CM))) |>
    ungroup()
}

v_aa <- aa_dat(valid_amino_acids[7])
table(v_aa$codon, v_aa$CM)

gln_state <- v_aa %>%
  filter(codon == "CAG") %>%
  select(organism_id, gln_state = CM)

v_aa <- v_aa %>%
  left_join(
    gln_state,
    by = "organism_id"
  )

ggplot(v_aa,
       aes(x = gln_state,
           y = RSCU)) +
  geom_boxplot() +
  facet_wrap(~ codon)

mod <- lm(
  RSCU ~ codon * gln_state +
    GC3 +
    phyPC1 + phyPC2 + phyPC3 + phyPC4,
  data = v_aa
)

anova(mod)

library(emmeans)

emm <- emmeans(mod, ~ codon | gln_state)

emm


emm_df <- as.data.frame(emm)

ggplot(emm_df,
       aes(codon, emmean,
           fill = gln_state)) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(ymin = lower.CL,
        ymax = upper.CL),
    position = position_dodge(0.9),
    width = 0.2
  )

v_aa %>%
  filter(codon == "CAG") %>%
  count(phylum, CM)

ggplot(
  v_aa %>%
    filter(codon == "CAG"),
  aes(phyPC1, fill = CM)
) +
  geom_density(alpha = 0.4)

ggplot(
  v_aa %>%
    filter(codon == "CAG"),
  aes(phyPC2, fill = CM)
) +
  geom_density(alpha = 0.4)
