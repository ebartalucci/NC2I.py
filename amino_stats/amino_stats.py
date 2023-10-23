# Plot statistics of amino acids in protein sequence
# Ettore Bartalucci, Aachen 20.10.23

import matplotlib.pyplot as plt
import numpy as np

def spacing_prot_seq(prot_seq, spaced_prot_seq):
    try:
        with open(prot_seq, 'r') as f:
            sequence = f.read()

        spaced = ' '.join(sequence)

        with open(spaced_prot_seq, 'w') as f:
            f.write(spaced)

        print('Protein sequence from Uniprot now contains spaces between each amino acid letter')

    except FileNotFoundError:
        print('Input file not found')
    except:
        print('An error occurred')

def count_amino_acids(prot_seq, count_file):
    try:
        with open(prot_seq, 'r') as f:
            sequence = f.read()

        amino_acid_count = {}
        for amino_acid in sequence:
            if amino_acid in amino_acid_count:
                amino_acid_count[amino_acid] += 1
            else:
                amino_acid_count[amino_acid] = 1

        with open(count_file, 'w') as f:
            for amino_acid, count in amino_acid_count.items():
                f.write(f"{amino_acid}: {count}\n")

        print('Amino acid counts written to amino_acid_count.txt')

    except FileNotFoundError:
        print('Input file not found')
    except:
        print('An error occurred')

def plot_amino_acid_statistics(count_file, plot_file):
    amino_acid_counts = {}
    total_count = 0

    with open(count_file, 'r') as f:
        for line in f:
            parts = line.strip().split(': ')
            if len(parts) == 2:
                amino_acid, count = parts
                amino_acid_counts[amino_acid] = int(count)
                total_count += int(count)

    amino_acids = list(amino_acid_counts.keys())
    counts = list(amino_acid_counts.values())

    percentages = [count / total_count * 100 for count in counts]

    plt.bar(amino_acids, percentages)
    plt.xlabel('Amino Acid')
    plt.ylabel('Percentage (%)')
    plt.ylim(0, 20)  # Set the y-axis limit to 0% to 18%
    plt.title('Amino Acid Occurrence Statistics')
    
    for i in range(len(amino_acids)):
        plt.text(amino_acids[i], percentages[i], f"{counts[i]} ({percentages[i]:.2f}%)", ha='center', va='bottom', rotation=90)
    
    plt.savefig(plot_file)
    plt.show()

protein_sequence = 'fus_fl_human.txt'
spaced_sequence = 'spaced_fus_fl_human.txt'
count_file = 'amino_acid_count.txt'
plot_file = 'amino_acid_statistics.png'

spacing_prot_seq(protein_sequence, spaced_sequence)
count_amino_acids(spaced_sequence, count_file)
plot_amino_acid_statistics(count_file, plot_file)
