import pandas as pd 
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description="Generate coverage heatmap.")
parser.add_argument("-i", "--input_dir", required=True)
parser.add_argument("-r", "--region", required=True)
parser.add_argument("-o", "--output_dir", required=True)
parser.add_argument("-n", "--num_samples", type=int, required=True)

args = vars(parser.parse_args())

input_dir = args["input_dir"]
region = args["region"].replace(":", "_").replace("-", "_")
output_dir = args["output_dir"]
num_samples = args["num_samples"]

barcode_ids = [f"{i:02}" for i in range(1, num_samples + 1)]
bin_size = 1000  # Bin size in bp

# group data by bin + mean depth in each bin 
def bin_coverage(df, bin_size):
    df = df.copy()
    start_pos = df["pos"].min() 
    df["bin"] = ((df["pos"] - start_pos) // bin_size) * bin_size + start_pos # makes bin alignment relative to the minimum position
    # integer division for position determination
    return df.groupby("bin")["depth"].mean().reset_index() 

# data matrix: rows = samples, columns = binned positions
coverage_matrix = []
bin_positions = None
sample_labels = []

for barcode in barcode_ids:
    file_path = os.path.join(input_dir, f"depth_barcode{barcode}.txt")
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        try:
            df = pd.read_csv(file_path, sep='\t', header=None, names=["chrom", "pos", "depth"])
            binned_df = bin_coverage(df, bin_size)
            if bin_positions is None:
                bin_positions = binned_df["bin"]
            coverage_matrix.append(binned_df["depth"].values)
            sample_labels.append(f"Sample {barcode}")
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
    else:
        print(f"File {file_path} does not exist or is empty.")

# Convert to 2D array (samples x bins)
coverage_array = np.array(coverage_matrix)

# plot heatmap
plt.figure(figsize=(15, 6)) 
sns.heatmap(
    coverage_array,
    xticklabels=False,  # Remove automatic tick labels
    yticklabels=sample_labels,
    cmap="BuPu",
    cbar_kws={'label': 'Coverage depth'}
)

interval = 10  # show every 6th bin
tick_positions = np.arange(0, len(bin_positions), interval)
tick_labels = bin_positions[::interval].astype(str)  # Get genomic positions at intervals

plt.xticks(tick_positions, tick_labels, rotation=90, ha="center")
plt.xlabel("Genomic position (1000 bp bins)") 
plt.title(f"Coverage heatmap: {args['region']}")
plt.tight_layout()
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, f"coverage_heatmap_{region}.png")
plt.savefig(output_file, dpi=300)
print(f"Saved heatmap to: {output_file}")
