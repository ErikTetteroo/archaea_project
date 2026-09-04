#!/bin/bash

# Output file
output="trnas_summary.csv"

# Write header
echo "organism,tRNA_type,anticodon,inf_score,note" > "$output"

# Loop through all trnas.out files
for file in *_trnas.out; do
    # Extract organism name (remove suffix)
    organism=$(basename "$file" _trnas.out)

    # Process file
    awk -v org="$organism" '
    BEGIN { OFS="," }
    # Skip header lines (first 3 lines and dashed separator)
    NR > 3 && $1 !~ /^-+/ {
        # Columns:
        # $5 = Type, $6 = Codon, $9 = Score, $10 = Note (may be empty)
        note = ($10 == "" ? "" : $10)
        print org, $5, $6, $9, note
    }
    ' "$file" >> "$output"

done