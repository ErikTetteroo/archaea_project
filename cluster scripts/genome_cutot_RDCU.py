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

codon_to_aa = {}
for aa, codons in codon_order:
    for c in codons:
        codon_to_aa[c] = aa

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
    codon_pair_counts = Counter()
    aa_counts = Counter()
    dipeptide_counts = Counter()

    # Collect codons from all genes
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        
        codons = [
            seq[i:i+3]
            for i in range(0, len(seq), 3)
            if len(seq[i:i+3]) == 3
        ]
        global_counts.update(codons)
        
        # Count amino acids
        aas = [
            codon_to_aa[c]
            for c in codons
            if c in codon_to_aa and codon_to_aa[c] != "TER"
        ]
        aa_counts.update(aas)

        # Count pairs
        for i in range(len(codons) - 1):
            c1, c2 = codons[i], codons[i+1]

            if (
                c1 in codon_to_aa and
                c2 in codon_to_aa and
                codon_to_aa[c1] != "TER" and
                codon_to_aa[c2] != "TER"
            ):
                codon_pair_counts[(c1, c2)] += 1

                a1 = codon_to_aa[c1]
                a2 = codon_to_aa[c2]
                dipeptide_counts[(a1, a2)] += 1

    total_codons = sum(global_counts.values())
    total_pairs = sum(codon_pair_counts.values())

    f_codon = {c: count / total_codons for c, count in global_counts.items()}
    total_aa = sum(aa_counts.values())
    f_aa = {a: count / total_aa for a, count in aa_counts.items()}
    
    E_aa = {}
    for a1 in f_aa:
        for a2 in f_aa:
            E_aa[(a1, a2)] = f_aa[a1] * f_aa[a2] * total_pairs

    norm = {}
    for (a1, a2), obs in dipeptide_counts.items():
        exp = E_aa.get((a1, a2), 0)
        if exp > 0:
            norm[(a1, a2)] = obs / exp
        else:
            norm[(a1, a2)] = 0

    RDCU = {}
    
    E_initial = {}
    for c1 in f_codon:
        for c2 in f_codon:
            E_initial[(c1, c2)] = f_codon[c1] * f_codon[c2] * total_pairs
    
    for (c1, c2), obs in codon_pair_counts.items():

        a1 = codon_to_aa[c1]
        a2 = codon_to_aa[c2]

        exp_init = E_initial.get((c1, c2), 0)
        coeff = norm.get((a1, a2), 0)

        exp_corrected = exp_init * coeff

        if exp_corrected > 0:
            RDCU[(c1, c2)] = obs / exp_corrected
        else:
            RDCU[(c1, c2)] = None

    # -----------------------------
    # Write RDCU output
    # -----------------------------
    rdcu_output_path = output_path.replace(".blk", "_RDCU.csv")

    with open(rdcu_output_path, "w") as out:

        # Header
        out.write(
        "codon_pair,c1,c2,RDCU,observed,expected,aa1,aa2,a_observed,a_expected\n"
        )

        for (c1, c2), obs in codon_pair_counts.items():

            a1 = codon_to_aa[c1]
            a2 = codon_to_aa[c2]

            # Expected codon pair (before normalization)
            exp_init = E_initial.get((c1, c2), 0)

            # Dipeptide observed/expected
            a_obs = dipeptide_counts.get((a1, a2), 0)
            a_exp = E_aa.get((a1, a2), 0)

            # RDCU
            rdcu_val = RDCU.get((c1, c2), None)

            # Handle None safely
            if rdcu_val is None:
                rdcu_str = ""
            else:
                rdcu_str = f"{rdcu_val:.6f}"

            # Write row
            out.write(
                f"{c1+c2},{c1},{c2},{rdcu_str},{obs},{exp_init:.6f},"
                f"{a1},{a2},{a_obs},{a_exp:.6f}\n"
            )

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