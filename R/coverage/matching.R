give_empty_codon_matrix <- function(
    codons,
    anticodons
) {
  coverage_mat <- matrix(
    0,
    nrow = length(anticodons),
    ncol = length(codons),
    dimnames = list(
      anticodon = anticodons,
      codon = codons
    )
  )
  coverage_mat
}

possible_codons <- function(anticodon, wobble_bases = NULL) {
  
  bases <- strsplit(anticodon, "")[[1]]
  
  # standard first two positions
  c1 <- EXACT_PAIRINGS[[bases[1]]]
  c2 <- EXACT_PAIRINGS[[bases[2]]]
  
  # wobble position
  c3 <- if (is.null(wobble_bases)) {
    EXACT_PAIRINGS[[bases[3]]]
  } else {
    wobble_bases
  }
  
  paste0(c1, c2, c3)
}

can_use_inosine <- function(
    anticodon, 
    codon_tbl
) {
  bases <- strsplit(anticodon, "")[[1]]
  
  # inosine derived from A-to-I editing
  if (bases[3] != "A") {
    return(FALSE)
  }
  
  aa <- codon_tbl %>%
    filter(anticodon == !!anticodon) %>%
    pull(aa) %>%
    first()
  
  codons <- possible_codons(
    anticodon,
    wobble_bases = INOSINE_PAIRINGS
  )
  
  all(
    codon_tbl %>%
      filter(codon_rna %in% codons) %>%
      pull(aa) == aa
  )
}

can_superwobble <- function(
    anticodon, 
    codon_tbl
)  {
    bases <- strsplit(anticodon, "")[[1]]
    
    # superwobble U at wobble position
    if (bases[3] != "U") {
      return(FALSE)
    }
    
    aa <- codon_tbl %>%
      filter(anticodon == !!anticodon) %>%
      pull(aa) %>%
      first()
    
    codons <- possible_codons(
      anticodon,
      wobble_bases = SUPERWOBBLE_PAIRINGS
    )
    
    all(
      codon_tbl %>%
        filter(codon_rna %in% codons) %>%
        pull(aa) == aa
    )
}

#function for checking what kind of pairing is possible for a specific codon and anticodon
check_pair <- function(
    codon,
    anticodon,
    inosine_allowed = FALSE,
    superwobble_allowed = FALSE,
    allow_inosine = FALSE,
    allow_superwobble = TRUE
) {
  
  codon_bases <- strsplit(codon, "")[[1]]
  anti_bases  <- strsplit(anticodon, "")[[1]]
  
  # positions 1 and 2 must pair perfectly
  for (i in 1:2) {
    
    if (EXACT_PAIRINGS[[anti_bases[i]]] != codon_bases[i]) {
      return(0)
    }
  }
  
  anti3  <- anti_bases[3]
  codon3 <- codon_bases[3]
  
  # exact pairing
  if (EXACT_PAIRINGS[[anti3]] == codon3) {
    return("M")
  }
  
  # GU wobble
  if (
    anti3 %in% names(GU_WOBBLE_PAIRINGS) &&
    GU_WOBBLE_PAIRINGS[[anti3]] == codon3
  ) {
    return("GU")
  }
  
  # inosine wobble
  if (
    allow_inosine &&
    inosine_allowed &&
    anti3 == "A" &&
    codon3 %in% INOSINE_PAIRINGS
  ) {
    return("I")
  }
  
  # superwobble
  if (
    allow_superwobble &&
    superwobble_allowed &&
    anti3 == "U" &&
    codon3 %in% SUPERWOBBLE_PAIRINGS
  ) {
    return("SU")
  }
  
  0
}