library(dplyr)

# Individual codon frequencies
codon_dat <- merged %>%
  filter(codon %in% c("CCG", "AAG")) %>%
  select(
    organism_id,
    codon,
    total,
    organism_total,
    RSCU,
    GC3,
    CM
  ) %>%
  mutate(
    observed_freq = total / organism_total
  ) %>%
  select(
    organism_id, codon, RSCU, GC3, CM, observed_freq
  ) %>%
  tidyr::pivot_wider(
    names_from = codon,
    values_from = c(RSCU, CM, observed_freq),
    names_sep = "_"
  )


# CCGAAG pair data
pair_dat <- pairs %>%
  filter(codon_pair == "CCGAAG") %>%
  mutate(
    pair_O_E = observed / expected,
    dipeptide_O_E = a_observed / a_expected,
    RDCU_check = pair_O_E / dipeptide_O_E
  ) %>%
  select(
    organism,
    codon_pair,
    GC3,
    CMP,
    observed,
    expected,
    a_observed,
    a_expected,
    RDCU,
    pair_O_E,
    dipeptide_O_E,
    RDCU_check
  )


# Combine
ccgaag_diag <- pair_dat %>%
  left_join(
    codon_dat,
    by = c("organism" = "organism_id")
  )


library(ggplot2)

# CCG individual preference
ggplot(
  ccgaag_diag,
  aes(GC3.x, RSCU_CCG, colour = CMP)
) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  labs(
    title = "CCG preference",
    y = "RSCU (CCG)",
    x = "GC3"
  )

ggplot(
  ccgaag_diag,
  aes(GC3.x, pair_O_E, colour = CMP)
) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  labs(
    title = "CCGAAG raw pair enrichment",
    y = "Observed / expected",
    x = "GC3"
  )

ggplot(
  ccgaag_diag,
  aes(GC3.x, RDCU, colour = CMP)
) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  labs(
    title = "CCGAAG RDCU",
    y = "RDCU",
    x = "GC3"
  )

ggplot(
  ccgaag_diag,
  aes(observed_freq_CCG, expected, colour = CMP)
) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  labs(
    title = "CCG frequency vs expected CCGAAG frequency",
    x = "Observed CCG frequency",
    y = "Expected CCGAAG frequency"
  )

ccgaag_diag %>%
  summarise(
    cor_CCG_RSCU_pair_OE =
      cor(RSCU_CCG, pair_O_E, use = "complete.obs"),
    
    cor_CCG_freq_pair_OE =
      cor(observed_freq_CCG, pair_O_E, use = "complete.obs"),
    
    cor_pair_OE_RDCU =
      cor(pair_O_E, RDCU, use = "complete.obs"),
    
    cor_dipeptide_RDCU =
      cor(dipeptide_O_E, RDCU, use = "complete.obs")
  )

ccgaag_diag %>%
  filter(RDCU > 3) %>%
  select(
    organism,
    GC3.x,
    CMP,
    observed,
    expected,
    pair_O_E,
    dipeptide_O_E,
    RDCU
  )

c(GCF_000621965.1, GCF_002287195.1, GCF_049554935.1)
unique(merged[merged$organism_id=="GCF_000621965.1",18:23])
unique(merged[merged$organism_id=="GCF_002287195.1",18:23])
unique(merged[merged$organism_id=="GCF_049554935.1",18:23])
