import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import sys

sns.set(style="whitegrid")
plt.rcParams.update({'figure.max_open_warning': 0})

bins_10 = np.arange(0, 1000, 10)
bins_50 = np.arange(1000, 3050, 50)  
bin = np.concatenate((bins_10, bins_50))

sample_num = sys.argv[1]

basecalled = f"ReadLengths/all_read_lengths_barcode{sample_num}.txt"
mapped = f"ReadLengths/all_mapped_read_lengths_barcode{sample_num}.txt"

try:
    basecalled_lengths = pd.read_csv(basecalled, header=None, names=["Length"])
    mapped_lengths = pd.read_csv(mapped, header=None, names=["Length"])
except FileNotFoundError as e:
    print(f"File not found: {e}")
    sys.exit(1)

basecalled_lengths = basecalled_lengths[basecalled_lengths["Length"] <= 3000]
mapped_lengths = mapped_lengths[mapped_lengths["Length"] <= 3000]

def plot_histogram(data, title, output_file):
    data = data.copy()
    data['Bin'] = pd.cut(data['Length'], bins=bin, right=False)

    binned = data['Bin'].value_counts().reset_index()
    binned.columns = ['Interval', 'Count']
    binned = binned.sort_values('Interval')

    plt.figure(figsize=(28, 10))
    sns.barplot(x='Interval', y='Count', data=binned, color="skyblue")

    plt.title(title)
    plt.xlabel("Read length intervals (bp)")
    plt.ylabel("Frequency")

    xticks = binned['Interval'].astype(str).tolist()
    step = 1  # Show all the labels on xaxis!
    plt.xticks(ticks=range(0, len(xticks), step), labels=xticks[::step], rotation=90, ha='center')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

os.makedirs("ReadLengths/Histogram", exist_ok=True)
plot_histogram(basecalled_lengths, f"Basecalled Read Lengths - Barcode {sample_num}", f"ReadLengths/Histogram/histogram_basecalled_barcode{sample_num}.png")
plot_histogram(mapped_lengths, f"Mapped Read Lengths - Barcode {sample_num}", f"ReadLengths/Histogram/histogram_mapped_barcode{sample_num}.png")
print(f"Histograms (basecalled & mapped) ok for barcode {sample_num}") 
