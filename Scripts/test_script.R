library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)

merged <- read_csv("Data/Merged_data_m.csv")
standard_set_coverage_table <- read_csv("Data/standard_set_coverage_table.csv")


merged_standard <- merged %>%
  
  # Remove old coverage columns
  select(-M, -GUw, -Iw, -CV) %>%
  
  # Add standard-set coverage values by codon
  left_join(
    standard_set_coverage_table %>%
      select(codon, M, GUw, CV),
    
    by = "codon"
  )




# Expand organism list-column
impossible_states <- possible_aa_delta_states %>%
  filter(impossible == 1) %>%
  unnest(organism_ids) %>%
  rename(
    organism_id = organism_ids
  )

merged_possible <- merged_standard %>%
  anti_join(
    impossible_states %>%
      select(organism_id, tRNA_type),
    
    by = c(
      "organism_id",
      "amino_acid" = "tRNA_type"
    )
  )

possible_states_only <- possible_aa_delta_states %>%
  filter(impossible == 0)

possible_expanded <- possible_states_only %>%
  unnest(organism_ids) %>%
  rename(
    organism_id = organism_ids
  )

parsed_deltas <- possible_expanded %>%
  
  separate_rows(aa_signature, sep = ";") %>%
  
  extract(
    aa_signature,
    into = c("codon", "delta_M", "delta_GUw", "delta_CV"),
    regex = "([AUGC]{3})\\(([-0-9.]+),([-0-9.]+),([-0-9.]+)\\)",
    convert = TRUE
  )

merged_updated <- merged_possible %>%
  
  left_join(
    parsed_deltas %>%
      select(
        organism_id,
        codon,
        delta_M,
        delta_GUw,
        delta_CV
      ),
    
    by = c("organism_id", "codon")
  ) %>%
  
  mutate(
    M   = M   + coalesce(delta_M, 0),
    GUw = GUw + coalesce(delta_GUw, 0),
    CV  = CV  + coalesce(delta_CV, 0)
  ) %>%
  
  select(
    -delta_M,
    -delta_GUw,
    -delta_CV
  )

merged_updated <- merged_updated %>%
  filter(
    amino_acid != "TER",
    !(amino_acid == "Met" & CV == 0)
  )

merged_updated$CV[merged_updated$CV==0] <- 0.4

arg <- merged_updated[merged_updated$amino_acid=="Arg",]
unique(arg$CV)


ala$codon <- factor(ala$codon)

ala_gcg$CV <- factor(ala_gcg$CV)

summary(
  lm(RSCU ~ CV + GC3, data = ala_gcg)
)


library(tidyr)
library(dplyr)

# Expand organism list-column
arg_states <- possible_states_only %>%
  
  filter(tRNA_type == "Arg") %>%
  
  unnest(organism_ids) %>%
  
  rename(
    organism_id = organism_ids
  ) %>%
  
  select(
    organism_id,
    aa_signature
  )

# Join onto Arg codon dataframe
arg <- arg %>%
  
  left_join(
    arg_states,
    by = "organism_id"
  ) %>%
  
  mutate(
    aa_signature = ifelse(
      is.na(aa_signature),
      "standard",
      aa_signature
    )
  )

arg$aa_signature <- factor(
  arg$aa_signature,
  levels = c(
    "standard",
    "CGG(-1,0,-1)",
    "AGG(-1,0,-1);CGG(-1,0,-1)",
    "AGG(-1,0,-1)",
    "CGG(-1,0,-1);CGU(1,0,1)"
  )
)

lm(RSCU ~ aa_signature * codon + GC3, data = arg)

arg_summary <- arg %>%
  group_by(aa_signature, codon) %>%
  summarise(
    mean_RSCU = mean(RSCU),
    mean_GC3 = mean(GC3),
    n = n()
  )

library(ggplot2)

ggplot(arg_summary,
       aes(codon, aa_signature, fill = mean_RSCU)) +
  geom_tile()
