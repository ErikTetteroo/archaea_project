#!/bin/bash

output="combined_codon_usage.csv"

echo "organism,codon,amino_acid,First,Second,Third,GC_third,organism_total,total,RSCU,GC3" > "$output"

for file in *.blk; do
    organism=$(basename "$file" .blk)

    # First pass: extract metadata
    GC3=$(grep "GC3s value in genome" "$file" | awk '{print $1}')
    total_codons=$(grep "total codons in genome" "$file" | awk '{print $1}')

    # Second pass: parse codon table
    awk -v org="$organism" -v GC3="$GC3" -v total_codons="$total_codons" '
    /^([A-Z][a-z]{2}|[A-Z]{3})[[:space:]]/ {
        aa = $1

        for (i = 2; i <= NF; i += 3) {
            codon = $i
            count = $(i+1)
            rscu  = $(i+2)

            if (codon == "" || count == "" || rscu == "")
                continue

            # Split codon into nucleotides
            first  = substr(codon, 1, 1)
            second = substr(codon, 2, 1)
            third  = substr(codon, 3, 1)

            # GC at third position
            if (third == "G" || third == "C")
                gc_third = 1
            else
                gc_third = 0

            print org "," codon "," aa "," first "," second "," third "," gc_third "," total_codons "," count "," rscu "," GC3
        }
    }
    ' "$file" >> "$output"

done

# -----------------------------
# Combine RDCU files
# -----------------------------
rdcu_output="combined_codon_pair_usage.csv"

echo "organism,codon_pair,c1,c2,RDCU,observed,expected,aa1,aa2,a_observed,a_expected" > "$rdcu_output"

for file in *_RDCU.csv; do
    organism=$(basename "$file" _RDCU.csv)

    awk -v org="$organism" 'NR > 1 {
        print org "," $0
    }' "$file" >> "$rdcu_output"

done