#!/usr/bin/env python3

import sys
import os
from collections import Counter
from Bio import SeqIO


# -----------------------------
# Codon order
# -----------------------------
codon_order = [
    ('Phe', ['TTT','TTC']),
    ('Leu', ['TTA','TTG','CTT','CTC','CTA','CTG']),
    ('Ile', ['ATT','ATC','ATA']),
    ('Met', ['ATG']),
    ('Val', ['GTT','GTC','GTA','GTG']),
    ('Ser', ['TCT','TCC','TCA','TCG','AGT','AGC']),
    ('Pro', ['CCT','CCC','CCA','CCG']),
    ('Thr', ['ACT','ACC','ACA','ACG']),
    ('Ala', ['GCT','GCC','GCA','GCG']),
    ('Tyr', ['TAT','TAC']),
    ('His', ['CAT','CAC']),
    ('Gln', ['CAA','CAG']),
    ('Asn', ['AAT','AAC']),
    ('Lys', ['AAA','AAG']),
    ('Asp', ['GAT','GAC']),
    ('Glu', ['GAA','GAG']),
    ('Cys', ['TGT','TGC']),
    ('Trp', ['TGG']),
    ('Arg', ['CGT','CGC','CGA','CGG','AGA','AGG']),
    ('Gly', ['GGT','GGC','GGA','GGG']),
    ('TER', ['TAA','TAG','TGA'])
]

# ==============================
# GC3 synonymous
# ==============================

def gc3_from_counter(codon_counts):
    exclude_codons = {"ATG", "TGG", "TAA", "TAG", "TGA"}

    gc = 0
    total = 0

    for codon, count in codon_counts.items():

        if codon in exclude_codons:
            continue

        total += count

        if codon[2] in {"G", "C"}:
            gc += count

    if total == 0:
        return None

    return gc / total
    
# -----------------------------
# Main function
# -----------------------------
def genome_cutot(fasta_file, output_path):

    global_counts = Counter()

    # Collect codons from all genes
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        codons = [
            seq[i:i+3]
            for i in range(0, len(seq), 3)
            if len(seq[i:i+3]) == 3
        ]
        global_counts.update(codons)
        
    gc3 = gc3_from_counter(global_counts)

    total_codons = sum(global_counts.values())

    with open(output_path, "w") as out:

        out.write("===== GENOME-WIDE CODON USAGE =====\n\n")

        for aa, codon_list in codon_order:

            total = sum(global_counts.get(c, 0) for c in codon_list)
            out.write(f"{aa:<4} ")

            for codon in codon_list:
                observed = global_counts.get(codon, 0)

                if total > 0:
                    rscu = observed / (total / len(codon_list))
                else:
                    rscu = 0.0

                codon_rna = codon.replace("T", "U")
                out.write(f"{codon_rna:>3} {observed:>6} {rscu:>6.2f} ")

            out.write("\n")

        out.write(f"\n{gc3:.2f} GC3s value in genome\n")
        out.write(f"\n{total_codons} total codons in genome\n")


# -----------------------------
# Command-line handling
# -----------------------------
if len(sys.argv) != 2:
    print("Usage: python genome_cutot.py <cds_fasta>")
    sys.exit(1)

fasta_file = sys.argv[1]

# Create output directory
output_dir = "codon_analysis_results_genomewise"
os.makedirs(output_dir, exist_ok=True)

# Create output filename
base = os.path.splitext(os.path.basename(fasta_file))[0]
output_path = os.path.join(output_dir, f"{base}.blk")

# Run analysis
genome_cutot(fasta_file, output_path)

print(f"{base} Finished")