import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import matplotlib.font_manager as fm
import os

# Load Roboto font from local file
script_dir = os.path.dirname(os.path.abspath(__file__))
roboto_path = os.path.join(script_dir, 'Roboto-Regular.ttf')
if os.path.exists(roboto_path):
    fm.fontManager.addfont(roboto_path)
    plt.rcParams['font.family'] = 'Roboto'
else:
    plt.rcParams['font.family'] = 'Arial'

# Read the Excel file
df = pd.read_excel(r'c:\Users\gliu6\OneDrive\Desktop\Engineered proteins\Attempt 3 LanM Li\LanM Li V3 results VBA enabled.xlsm', sheet_name='mpnn_results')

# Get sequences
sequences = df['seq'].tolist()

# Find variable positions
seq_len = len(sequences[0])
position_variation = []
for pos in range(seq_len):
    residues = [seq[pos] for seq in sequences if len(seq) > pos]
    unique = len(set(residues))
    position_variation.append((pos, unique))

variable_positions = [p for p, v in position_variation if v > 1]

# Group consecutive variable positions
groups = []
current_group = [variable_positions[0]]
for pos in variable_positions[1:]:
    if pos - current_group[-1] <= 3:
        current_group.append(pos)
    else:
        groups.append(current_group)
        current_group = [pos]
groups.append(current_group)

# Define amino acids
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def calculate_prevalence(sequences, start_pos, end_pos):
    num_positions = end_pos - start_pos
    prevalence = np.zeros((len(amino_acids), num_positions))
    for pos_idx, pos in enumerate(range(start_pos, end_pos)):
        residues = [seq[pos] for seq in sequences if len(seq) > pos]
        counts = Counter(residues)
        total = len(residues)
        for aa_idx, aa in enumerate(amino_acids):
            prevalence[aa_idx, pos_idx] = (counts.get(aa, 0) / total) * 100
    return prevalence

# Get loop regions
loop1_start, loop1_end = min(groups[0]), max(groups[0]) + 1
loop2_start, loop2_end = min(groups[1]), max(groups[1]) + 1

print(f"Loop 1: positions {loop1_start}-{loop1_end-1}")
print(f"Loop 2: positions {loop2_start}-{loop2_end-1}")

# Calculate prevalence
prevalence1 = calculate_prevalence(sequences, loop1_start, loop1_end)
prevalence2 = calculate_prevalence(sequences, loop2_start, loop2_end)

# Combine the two heatmaps side by side with a gap
gap = np.full((len(amino_acids), 1), np.nan)  # Gap column
combined = np.hstack([prevalence1, gap, prevalence2])

# Create compact figure with shared Y-axis and single colorbar
fig, ax = plt.subplots(figsize=(14, 8), dpi=300)

# Custom colormap
cmap = sns.color_palette("YlOrBr", as_cmap=True)
cmap.set_bad(color='white')  # White for the gap

# Create x-axis labels
n_pos1 = prevalence1.shape[1]
n_pos2 = prevalence2.shape[1]
x_labels = [str(i) for i in range(1, n_pos1 + 1)] + [''] + [str(i) for i in range(1, n_pos2 + 1)]

# Plot combined heatmap
hm = sns.heatmap(combined, ax=ax, cmap=cmap, vmin=0, vmax=100,
                 xticklabels=x_labels,
                 yticklabels=amino_acids,
                 annot=True, fmt='.0f', annot_kws={'size': 7},
                 cbar_kws={'label': 'Prevalence (%)', 'shrink': 0.8})

# Add section labels with full descriptions
ax.text(n_pos1 / 2, -0.8, 'EF Hand 2 (Loop 1)\nAmino Acid Prevalence (%)', ha='center', fontsize=11, fontweight='bold')
ax.text(n_pos1 + 1 + n_pos2 / 2, -0.8, 'EF Hand 3 (Loop 2)\nAmino Acid Prevalence (%)', ha='center', fontsize=11, fontweight='bold')

ax.set_xlabel('Position', fontsize=12)
ax.set_ylabel('Amino Acid', fontsize=12)

plt.tight_layout()

# Save
output_path = r'c:\Users\gliu6\OneDrive\Desktop\yayasoftware\ee360t\pst7\peptide_heatmap_compact.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nCompact heatmap saved to: {output_path}")

pdf_path = r'c:\Users\gliu6\OneDrive\Desktop\yayasoftware\ee360t\pst7\peptide_heatmap_compact.pdf'
plt.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"PDF version saved to: {pdf_path}")

plt.close()
print("\nDone!")
