#!/bin/bash

DATASETS="/home2/s4188829/Codon_project/programmes/datasets"
OUTDIR="/scratch/s4188829/genomes"

mkdir -p $OUTDIR

while read acc; do
    echo "Downloading $acc"
    
    $DATASETS download genome accession $acc \
        --include cds \
        --filename ${OUTDIR}/${acc}.zip

done < accessions.txt