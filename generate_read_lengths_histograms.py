import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# Plot settings
sns.set(style="whitegrid")
plt.rcParams.update({'figure.max_open_warning': 0})

# Generate custom bin edges
bins_10 = np.arange(0, 1000, 10)
bins_50 = np.arange(1000, 3050, 50)  # 2000 inclusive
bin_edges = np.concatenate((bins_10, bins_50))

# Loop through barcodes 01 to 12
for i in range(1, 13):
    barcode = f"{i:02d}"

    basecalled_file = f"all_read_lengths_barcode{barcode}.txt"
    mapped_file = f"all_mapped_read_lengths_barcode{barcode}.txt"

    try:
        basecalled_lengths = pd.read_csv(basecalled_file, header=None, names=["Length"])
        mapped_lengths = pd.read_csv(mapped_file, header=None, names=["Length"])
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        continue

    basecalled_lengths = basecalled_lengths[basecalled_lengths["Length"] <= 3000]
    mapped_lengths = mapped_lengths[mapped_lengths["Length"] <= 3000]

    def plot_histogram(data, title, output_file):
        data = data.copy()
        data['Bin'] = pd.cut(data['Length'], bins=bin_edges, right=False)

        binned = data['Bin'].value_counts().reset_index()
        binned.columns = ['Interval', 'Count']
        binned = binned.sort_values('Interval')

        plt.figure(figsize=(28, 10))
        sns.barplot(x='Interval', y='Count', data=binned, color="steelblue")

        plt.title(title)
        plt.xlabel("Read Length Intervals (bp)")
        plt.ylabel("Frequency")

        
        xticks = binned['Interval'].astype(str).tolist()
        step = 1  # Show every label
        plt.xticks(ticks=range(0, len(xticks), step), labels=xticks[::step], rotation=90, ha='center')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()

    plot_histogram(basecalled_lengths, f"Basecalled Read Lengths - Barcode {barcode}",
                   f"histogram_basecalled_barcode{barcode}.png")
    plot_histogram(mapped_lengths, f"Mapped Read Lengths - Barcode {barcode}",
                   f"histogram_mapped_barcode{barcode}.png")

print("Histograms generated for all barcodes.")
