# replace T with U
dna_to_rna <- function(x) gsub("T", "U", x)

# replace U with T
rna_to_dna <- function(x) gsub("U", "T", x)

# complementary sequence (not reversed)
comp <- function(seq) {
  chartr("AUCG", "UAGC", seq)
}

# reverse a sequence
reverse_seq <- function(seq) {
  sapply(strsplit(seq, ""), function(bases) {
    paste(rev(bases), collapse = "")
  })
}

# reverse + RNA->DNA formatting
reverse_complement_format <- function(x) {
  rna_to_dna(reverse_seq(x))
}

trnascan_to_anticodon <- function(x) {
  dna_to_rna(reverse_seq(x))
}