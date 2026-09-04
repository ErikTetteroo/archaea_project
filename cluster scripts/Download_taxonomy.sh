#!/bin/bash

DATASETS="/home2/s4188829/Codon_project/programmes/datasets"

echo -e "accession\torganism\ttaxid\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" > lineages.tsv

split -l 20 accessions.txt batch_

total=$(ls batch_* | wc -l)
current=0

for file in batch_*; do
  current=$((current + 1))
  echo "Batch $current / $total"

  # STEP 1: get accession, organism, taxid
  $DATASETS summary genome accession $(tr '\n' ' ' < "$file") \
  | jq -r '.reports[] | [.accession, .organism.organism_name, .organism.tax_id] | @tsv' \
  | while IFS=$'\t' read acc org taxid; do

      # STEP 2: get taxonomy for each taxid
      taxonomy=$($DATASETS summary taxonomy taxon "$taxid" \
        | jq -r '
          .reports[0].taxonomy.classification |
          [
            .domain.name // "NA",
            .phylum.name // "NA",
            .class.name // "NA",
            .order.name // "NA",
            .family.name // "NA",
            .genus.name // "NA",
            .species.name // "NA"
          ] | @tsv
        ')

      echo -e "$acc\t$org\t$taxid\t$taxonomy"

  done >> lineages.tsv

  sleep 1

done

rm batch_*