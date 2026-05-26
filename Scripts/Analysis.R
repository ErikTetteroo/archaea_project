#read in data
trnadata <- read.table("Data\\Raw_data\\trnas_summary.csv",sep = ",",header = T)


#filter to standard aa
aa <- c("Ala","Arg","Asn","Asp","Cys",
        "Gln","Glu","Gly","His","Ile",
        "Leu","Lys","Met","Phe","Pro",
        "Ser","Thr","Trp","Tyr","Val")

trnas_2 <- trnadata[trnadata$tRNA_type %in% aa,]

#filter out pseudogenes
trnas_2 <- trnadata[trnadata$note!="pseudo",]


n_occurtrna <- data.frame(table(trnas_2$organism))[]

n_occurinfo <- data.frame()
for (i in 1:length(unique(trnas_2$organism))) {
  n_occurinfo <- rbind(n_occurinfo,data.frame(organism = unique(trnas_2$organism)[i],
                                              table(trnas_2$anticodon[trnas_2$organism==unique(trnas_2$organism)[i]]))[])
}

n_occurorganism <- data.frame(table(n_occurinfo$organism))[]
n_occuranticodon <- data.frame(table(n_occurinfo$Var1))[]

#get some diagnostics
boxplot(n_occurorganism$Freq)
hist(n_occurinfo$Freq[n_occurinfo$Freq>1])

length(n_occurorganism$Freq[n_occurorganism$Freq<42])

#extract data on aa sets############################################################
library(dplyr)

result <- trnas_2 %>%
  group_by(organism) %>%
  summarise(anticodon = paste(sort(unique(anticodon)), collapse = ", ")) 
length(unique(result$anticodon[114]))

n_occursets <- data.frame(table(result$anticodon))[]

#order by set occurence
setssorted <- n_occursets[order(-n_occursets$Freq), ]
setssorted$setnr <- seq_len(nrow(setssorted))

result$setnr <- 0
for (i in 1:length(unique(result$anticodon))) {
  result$setnr[result$anticodon==unique(result$anticodon)[i]] <- setssorted$setnr[setssorted$Var1==unique(result$anticodon)[i]]
}

trnas_2$setnr <- 0
for (i in 1:length(unique(result$organism))) {
  trnas_2$setnr[trnas_2$organism==unique(result$organism)[i]] <- result$setnr[result$organism==unique(result$organism)[i]]
}



#create codon table
aatable <- data.frame(Firstcodon = c(rep("T",16),rep("C",16),rep("A",16),rep("G",16)),
                      Secondcodon = rep(c("T", "C", "A", "G"), each = 4),
                      Thirdcodon = rep(c("T", "C", "A", "G"), times = 8))
codons <- paste0(aatable$Firstcodon,aatable$Secondcodon,aatable$Thirdcodon)


#create heatmap
library(tidyr)
trnas_2o <- trnas_2[order(trnas_2$setnr),]

trnas_2o <- trnas_2o %>%
  group_by(setnr, organism) %>%
  summarise(n_dupes = n(), .groups = "drop") %>%   # count duplicates per organism in each set
  group_by(setnr) %>%
  arrange(setnr, desc(n_dupes)) %>%              # sort organisms by duplicates within each set
  mutate(setnr_rank = row_number()) %>%          # assign rank
  left_join(trnas_2o, by = c("setnr", "organism")) %>%
  ungroup()

trnas_2o <- trnas_2o %>%
  arrange(setnr, desc(n_dupes))

organism_order <- trnas_2o %>%
  distinct(organism, setnr)

df_counts <- expand.grid(
  organism = unique(trnas_2o$organism),
  anticodon = codons,
  stringsAsFactors = FALSE
) %>%
  left_join(
    trnas_2o %>% count(organism, anticodon),
    by = c("organism", "anticodon")
  ) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  left_join(trnas_2o %>% distinct(organism, setnr), by = "organism") %>%
  pivot_wider(names_from = anticodon, values_from = n) %>%
  arrange(setnr)

df_counts <- df_counts[,-2]

df_long <- df_counts %>%
  pivot_longer(
    cols = -organism,
    names_to = "codon",
    values_to = "amount"
  )

df_long$organism <- factor(df_long$organism, levels = rev(df_counts$organism))
df_long$codon <- factor(df_long$codon, levels = codons)

df_long <- df_long %>%
  mutate(organism_pos = as.numeric(factor(organism, levels = rev(unique(organism))))) # top-to-bottom
library(ggplot2)

ggplot(df_long, aes(x = codon, y = organism_pos, fill = amount)) +
  geom_tile() +
  scale_y_continuous(
    name = NULL,
    breaks = c(0, max(df_long$organism_pos)*0.25, max(df_long$organism_pos)*0.5,
               max(df_long$organism_pos)*0.75, max(df_long$organism_pos)),
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  scale_fill_gradientn(
    colors = c("white", "steelblue", "yellow"),
    values = scales::rescale(c(0, 1, 12))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_line() # show tick marks for percentages
  ) +
  ggtitle("tRNA pools of 764 Archaea genomes")
unique(df_long$organism)
