# pairing rules
EXACT_PAIRINGS <- c(
  A = "U",
  U = "A",
  G = "C",
  C = "G"
)

GU_WOBBLE_PAIRINGS <- c(
  G = "U",
  U = "G"
)

INOSINE_PAIRINGS <- c("U", "C", "A")

SUPERWOBBLE_PAIRINGS <- c("U", "C", "G")

# codon table
give_codon_table <- function()
{
  codon_tbl <- tribble(
    ~aa,   ~codons,
    "Phe", c("TTT","TTC"),
    "Leu", c("TTA","TTG","CTT","CTC","CTA","CTG"),
    "Ile", c("ATT","ATC","ATA"),
    "Met", c("ATG"),
    "Val", c("GTT","GTC","GTA","GTG"),
    "Ser", c("TCT","TCC","TCA","TCG","AGT","AGC"),
    "Pro", c("CCT","CCC","CCA","CCG"),
    "Thr", c("ACT","ACC","ACA","ACG"),
    "Ala", c("GCT","GCC","GCA","GCG"),
    "Tyr", c("TAT","TAC"),
    "His", c("CAT","CAC"),
    "Gln", c("CAA","CAG"),
    "Asn", c("AAT","AAC"),
    "Lys", c("AAA","AAG"),
    "Asp", c("GAT","GAC"),
    "Glu", c("GAA","GAG"),
    "Cys", c("TGT","TGC"),
    "Trp", c("TGG"),
    "Arg", c("CGT","CGC","CGA","CGG","AGA","AGG"),
    "Gly", c("GGT","GGC","GGA","GGG"),
    "TER", c("TAA","TAG","TGA")
  )

  codon_tbl <- codon_tbl %>%
  unnest(codons) %>%
  rename(codon_dna = codons) %>%
  mutate(
    codon_rna = dna_to_rna(codon_dna),
    anticodon = comp(codon_rna)
  )
  codon_tbl
}