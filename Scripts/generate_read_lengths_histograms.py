import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import sys

sns.set(style="whitegrid")
plt.rcParams.update({'figure.max_open_warning': 0})

sample_num = sys.argv[1]
workdir = sys.argv[2]
basecalled_file = os.path.join(workdir, f"ReadLengths/all_read_lengths_barcode{sample_num}.txt")
mapped_file = os.path.join(workdir, f"ReadLengths/all_mapped_read_lengths_barcode{sample_num}.txt")

try:
    basecalled_lengths = pd.read_csv(basecalled_file, header=None, names=["Length"])
    mapped_lengths = pd.read_csv(mapped_file, header=None, names=["Length"])
except FileNotFoundError as e:
    print(f"File not found: {e}")
    sys.exit(1)
basecalled_lengths = basecalled_lengths[basecalled_lengths["Length"] <= 3000]
mapped_lengths = mapped_lengths[mapped_lengths["Length"] <= 3000]
bins_10 = np.arange(0, 600, 10)
bins_50 = np.arange(600, 1005, 50) 
bin_edges = np.concatenate((bins_10, bins_50))

def plot_histogram(data, title, output_file):
    data = data[data['Length'] <= 1000]
    data['Bin'] = pd.cut(data['Length'], bins=bin_edges, right=False)
    binned = data['Bin'].value_counts().sort_index().reset_index()
    binned.columns = ['Interval', 'Count']
    binned['BinEnd'] = binned['Interval'].apply(lambda x: int(x.right))
    plt.figure(figsize=(18, 6))
    sns.barplot(x='BinEnd', y='Count', data=binned, color="blue")
    xticks = binned['BinEnd']
    xtick_labels = [str(x) if (x <= 600 and x % 50 == 0) or (x > 600 and x % 200 == 0) else "" for x in xticks]
    plt.xticks(ticks=range(len(xticks)), labels=xtick_labels, rotation=0, ha='right', fontsize=10)
    plt.title(title)
    plt.xlabel("Read length (bp)")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

output_dir = os.path.join(workdir, "ReadLengths/Histogram")
os.makedirs(output_dir, exist_ok=True)

plot_histogram(basecalled_lengths, f"Basecalled read lengths - Sample {sample_num}", os.path.join(output_dir, f"histogram_basecalled_barcode{sample_num}.png"))

plot_histogram(mapped_lengths, f"Mapped read lengths - Sample {sample_num}", os.path.join(output_dir, f"histogram_mapped_barcode{sample_num}.png"))

print(f"Histograms generated successfully for barcode {sample_num}")
