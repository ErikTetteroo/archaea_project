#load packages
library(dplyr)
library(stringr)
library(ggplot2)

#Data prep#####################################################################
#read in codon usage, trna data, and conversion table
Coverage_table <- read.table("Data\\Coverage_table.csv", sep = ',',header = T)[,-1]
Codon_usage <- read.table("Data\\Raw_data\\combined_codon_usage.csv", sep = ',', header = T)
tRNA_dat <- read.table("Data\\Raw_data\\trnas_summary.csv", sep = ',',header = T)
Taxa <- read.table("Data\\Raw_data\\lineages.tsv",sep = "\t", header = T)


#Create consistent organism id
trna <- tRNA_dat %>%
  mutate(organism_id = str_extract(organism, "GCF_[0-9]+\\.[0-9]+"))

#Filter trna pseudo genes & low score genes unrecognized anticodons
trna_filtered <- trna %>%
  filter(
    inf_score > 30,
    is.na(note) | !str_detect(note, "pseudo"),
    !str_detect(anticodon, "N")   
  )

#Merge trna & codon usage by id
codon_trna_pairs <- Codon_usage %>%
  rename(organism_id = organism) %>%
  inner_join(trna_filtered, by = "organism_id")

#Merge with coverage table by codon & anticodon
codon_trna_pairs <- codon_trna_pairs %>%
  inner_join(Coverage_table, by = c("codon", "anticodon"))

#Count pairing coverage each codon
pairing_counts <- codon_trna_pairs %>%
  filter(value != "0") %>%
  mutate(pair_type = case_when(
    value == "M"  ~ "M",
    value == "GU" ~ "GUw",
    value == "I"  ~ "Iw"
  )) %>%
  group_by(organism_id, codon, pair_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = pair_type,
    values_from = count,
    values_fill = 0
  )

#add pairing counts to final dataset
final_data <- Codon_usage %>%
  rename(organism_id = organism) %>%
  left_join(pairing_counts, by = c("organism_id", "codon")) %>%
  mutate(
    M   = coalesce(M, 0),
    GUw = coalesce(GUw, 0),
    Iw  = coalesce(Iw, 0)
  )

#combine pairing counts in Coverage value
GUweight <- 0.5
Iweight <- 0.5
final_data <- final_data %>%
  mutate(
    CV = M + GUweight * GUw + Iweight * Iw
  )

#add taxonomy
colnames(Taxa)[1] <- "organism_id"
Taxac <- Taxa[,-3:-4]

final_data_t <- left_join(final_data, Taxac, by = "organism_id")

write.csv(final_data_t, "Data\\Merged_data_m.csv")

###########################################################################

#Extracting special cases
selenocysteine_cases <- trna_filtered$organism_id[trna_filtered$tRNA_type=="SeC"]
sup_cases <- trna_filtered$organism_id[trna_filtered$tRNA_type=="Sup"]
inosine_cases <- final_data_t$organism_id[final_data_t$Iw>0]
uau <- final_data_t[final_data_t$codon=='AUA',]
uau_cases <- uau[uau$CV>0,]


###########################################################################


ggplot(final_data_t, aes(x = CV, y = RSCU)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~ amino_acid, scales = "free") +
  theme_minimal() +
  labs(
    title = "RSCU vs Coverage per amino acid"
  )

cor_summary <- final_data %>%
  group_by(amino_acid) %>%
  summarise(cor = cor(CV, RSCU, use = "complete.obs"))

cor_summary