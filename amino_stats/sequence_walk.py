import matplotlib.pyplot as plt

# Function to calculate amino acid occurrence statistics
def calculate_amino_acid_statistics(sequence):
    amino_acids = "ARNDCEQGHILKMFPSTWYV"  # List of 20 standard amino acids
    statistics = {aa: sequence.count(aa) for aa in amino_acids}
    return statistics

# Read the protein sequence from a text file
with open("fus_fl_human.txt", "r") as file:
    protein_sequence = file.read()

# Define domain positions
lc_start, lc_end = 0, 164
rgg1_start, rgg1_end = 165, 268  # Modified RGG1 start and end positions
rrm_start, rrm_end = 280, 376  # Modified RRM start and end positions
rgg2_start, rgg2_end = 377, 417  # Modified RGG2 start and end positions
rgg3_start, rgg3_end = 454, 500  # Modified RGG3 start and end positions
znf_start, znf_end = 418, 453  # Modified ZNF start and end positions

# Extract sequences for each domain
lc_sequence = protein_sequence[lc_start:lc_end]  # Extract sequence from lc_start to lc_end
rgg1_sequence = protein_sequence[rgg1_start:rgg1_end]  # Extract sequence from rgg1_start to rgg1_end
rrm_sequence = protein_sequence[rrm_start:rrm_end]  # Extract sequence from rrm_start to rrm_end
rgg2_sequence = protein_sequence[rgg2_start:rgg2_end]  # Extract sequence from rgg2_start to rgg2_end
rgg3_sequence = protein_sequence[rgg3_start:rgg3_end]  # Extract sequence from rgg3_start to rgg3_end
znf_sequence = protein_sequence[znf_start:znf_end]  # Extract sequence from znf_start to znf_end

# Calculate amino acid statistics for each domain
lc_stats = calculate_amino_acid_statistics(lc_sequence)
rgg1_stats = calculate_amino_acid_statistics(rgg1_sequence)
rrm_stats = calculate_amino_acid_statistics(rrm_sequence)
rgg2_stats = calculate_amino_acid_statistics(rgg2_sequence)
rgg3_stats = calculate_amino_acid_statistics(rgg3_sequence)
znf_stats = calculate_amino_acid_statistics(znf_sequence)

# Create subplots for LC, RGG1, RRM, RGG2, RGG3, and ZNF domains
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# LC subplot
axes[0, 0].bar(lc_stats.keys(), lc_stats.values(), label="LC")
axes[0, 0].set_title("Amino Acid Distribution in LC Domain")

# RGG1 subplot
axes[0, 1].bar(rgg1_stats.keys(), rgg1_stats.values(), label="RGG1")
axes[0, 1].set_title("Amino Acid Distribution in RGG1 Domain")

# RRM subplot
axes[0, 2].bar(rrm_stats.keys(), rrm_stats.values(), label="RRM")
axes[0, 2].set_title("Amino Acid Distribution in RRM Domain")

# RGG2 subplot
axes[1, 0].bar(rgg2_stats.keys(), rgg2_stats.values(), label="RGG2")
axes[1, 0].set_title("Amino Acid Distribution in RGG2 Domain")

# RGG3 subplot
axes[1, 1].bar(rgg3_stats.keys(), rgg3_stats.values(), label="RGG3")
axes[1, 1].set_title("Amino Acid Distribution in RGG3 Domain")

# ZNF subplot
axes[1, 2].bar(znf_stats.keys(), znf_stats.values(), label="ZNF")
axes[1, 2].set_title("Amino Acid Distribution in ZNF Domain")

# Set common y-axis label
for ax_row in axes:
    for ax in ax_row:
        ax.set_xlabel("Amino Acids")
        ax.set_ylabel("Count")

# Create a separate figure for the OUT domain
out_sequence = protein_sequence[lc_end:rgg1_start] + protein_sequence[rrm_end:rgg2_start] + protein_sequence[rgg3_end:]
out_stats = calculate_amino_acid_statistics(out_sequence)
fig_out, ax_out = plt.subplots(figsize=(6, 5))
ax_out.bar(out_stats.keys(), out_stats.values(), label="OUT")
ax_out.set_title("Amino Acid Distribution outside specified Domains")
ax_out.set_xlabel("Amino Acids")
ax_out.set_ylabel("Count")

# Calculate the total count of all amino acids
total_count = sum(lc_stats.values()) + sum(rgg1_stats.values()) + sum(rrm_stats.values()) + sum(rgg2_stats.values()) + \
    sum(rgg3_stats.values()) + sum(znf_stats.values()) + sum(out_stats.values())

# Create a final text annotation to display the total count
fig_total = plt.figure(figsize=(10, 6))
ax_total = fig_total.add_subplot(111)
ax_total.text(0.5, 0.5, f"Total Amino Acid Count: {total_count}", fontsize=16, ha="center")
ax_total.axis("off")

plt.tight_layout()
plt.show()
