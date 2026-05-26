#load packages
library(reshape2)
library(ggplot2)

#function for checking what kind of pairing is possible for a specific codon and anticodon
check_pair <- function(codon, anticodon, inosine_allowed) {
  codon_bases <- strsplit(codon, "")[[1]]
  anticodon_bases <- strsplit(anticodon, "")[[1]]
  
  # Positions 1 & 2 must match perfectly
  for (i in 1:2) {
    if (perfect_pairs[[anticodon_bases[i]]] != codon_bases[i]) {
      return(0)
    }
  }
  
  anti3 <- anticodon_bases[3]
  codon3 <- codon_bases[3]
  
  # Perfect
  if (!is.null(perfect_pairs[[anti3]]) && perfect_pairs[[anti3]] == codon3) {
    return(1)
  }
  
  # Standard wobble
  if (!is.null(wobble_pairs[[anti3]]) && codon3 %in% wobble_pairs[[anti3]]) {
    return(0.5)
  }
  
  # Inosine (ONLY if allowed for this anticodon)
  if (anti3 == "A" && inosine_allowed[anticodon] && codon3 %in% inosine_pairs) {
    return(0.5)
  }
  
  return(0)
}

#function for checking which A-ending codons cannot cause inosine wobbles
can_use_inosine <- function(anticodon, codon2aa) {
  
  bases <- strsplit(anticodon, "")[[1]]
  anti3 <- bases[3]
  
  # Only relevant if wobble base is A
  if (anti3 != "A") {
    return(FALSE)
  }
  
  # Get amino acid of this anticodon
  aa <- anticodon2aa[anticodon]
  
  # Generate all possible codons this anticodon could bind with inosine
  possible_codons <- c()
  
  for (b3 in c("U","C","A")) {  # inosine pairing
    codon <- paste0(
      perfect_pairs[[bases[1]]],
      perfect_pairs[[bases[2]]],
      b3
    )
    possible_codons <- c(possible_codons, codon)
  }
  
  # Check if all map to same amino acid
  all(codon2aa[possible_codons] == aa)
}

#function replacing T's with U's in codons
dna_to_rna <- function(x) gsub("T", "U", x)

#function for creating complementary codon/anticodon (does not reverse)
comp <- function(seq) {
  chartr("AUCG", "UAGC", paste(strsplit(seq, "")[[1]], collapse=""))
}

#function for returning anticodons to trnascan-se format
#(so reversed + replacing U's with T's)
reverse_complement_format <- function(x) {
  sapply(x, function(seq) {
    bases <- strsplit(seq, "")[[1]]
    reversed <- paste(rev(bases), collapse = "")
    gsub("U", "T", reversed)
  })
}

#############################################################################

#relate codons to amino acids
codon_order <- list(
  Phe = c("TTT","TTC"),
  Leu = c("TTA","TTG","CTT","CTC","CTA","CTG"),
  Ile = c("ATT","ATC","ATA"),
  Met = c("ATG"),
  Val = c("GTT","GTC","GTA","GTG"),
  Ser = c("TCT","TCC","TCA","TCG","AGT","AGC"),
  Pro = c("CCT","CCC","CCA","CCG"),
  Thr = c("ACT","ACC","ACA","ACG"),
  Ala = c("GCT","GCC","GCA","GCG"),
  Tyr = c("TAT","TAC"),
  His = c("CAT","CAC"),
  Gln = c("CAA","CAG"),
  Asn = c("AAT","AAC"),
  Lys = c("AAA","AAG"),
  Asp = c("GAT","GAC"),
  Glu = c("GAA","GAG"),
  Cys = c("TGT","TGC"),
  Trp = c("TGG"),
  Arg = c("CGT","CGC","CGA","CGG","AGA","AGG"),
  Gly = c("GGT","GGC","GGA","GGG"),
  TER = c("TAA","TAG","TGA")
)

codon2aa <- setNames(
  rep(names(codon_order), lengths(codon_order)),
  dna_to_rna(unlist(codon_order))
)

#relate anticodons to amino acids
anticodon2aa <- setNames(
  codon2aa,
  sapply(names(codon2aa), comp)
)

#check what they look like
codon2aa
anticodon2aa

#create pairing matrix
aatable <- data.frame(Firstcodon = c(rep("U",16),rep("C",16),rep("A",16),rep("G",16)),
                      Secondcodon = rep(c("U", "C", "A", "G"), each = 4),
                      Thirdcodon = rep(c("U", "C", "A", "G"), times = 8))
codons <- paste0(aatable$Firstcodon,aatable$Secondcodon,aatable$Thirdcodon)

coveragetable <- data.frame(codon = codons)
for (i in 1:length(codons)) {
  coveragetable <- cbind(coveragetable,0)
}
colnames(coveragetable) <- c('codon',codons)

#pairing logic
perfect_pairs <- list(
  A = "U",
  U = "A",
  G = "C",
  C = "G"
)

wobble_pairs <- list(
  G = "U",
  U = "G"
)

inosine_pairs <- "N" #c("U", "C", "A") turned off for now as inosine is not present in Archaea

inosine_allowed <- sapply(colnames(coveragetable)[-1], can_use_inosine, codon2aa = codon2aa)


#apply possible pairings to matrix
for (i in 1:nrow(coveragetable)) {
  codon <- coveragetable$codon[i]
  
  for (j in 2:ncol(coveragetable)) {
    anticodon <- colnames(coveragetable)[j]
    
    # Only allow same amino acid
    if (codon2aa[codon] != anticodon2aa[anticodon]) {
      next
    }
    
    coveragetable[i, j] <- check_pair(codon, anticodon, inosine_allowed)
  }
}

#convert to long format
coveragetable_mat <- coveragetable
rownames(coveragetable_mat) <- coveragetable_mat$codon
coveragetable_mat$codon <- NULL
df_long <- melt(cbind(codon = rownames(coveragetable_mat), coveragetable_mat),
                id.vars = "codon",
                variable.name = "anticodon",
                value.name = "value")

#Output
df_out <- df_long
df_out$anticodon <- as.character(df_out$anticodon)
df_out$anticodon <- reverse_complement_format(df_out$anticodon)

#write.csv(df_out, file = "Data\\Coverage_table.csv")


#Visualization plot
ggplot(df_long, aes(x = anticodon, y = codon, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", "lightblue", "blue"),
                       values = c(0, 0.5, 1),
                       limits = c(0,1.5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6))

