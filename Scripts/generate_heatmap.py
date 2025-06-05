import pandas as pd 
import numpy as np
import os, sys
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

ap = argparse.ArgumentParser(description="coverage heatmap")
ap.add_argument("-i", "--input_dir", required=True)
ap.add_argument("-r", "--region", required=True)
ap.add_argument("-o", "--output_dir", required=True)
ap.add_argument("-n", "--num_samples", type=int, required=True)
args = vars(ap.parse_args())

input = args["input_dir"]
region = args["region"].replace(":", "_").replace("-", "_") 
output = args["output_dir"]
num_samples = args["num_samples"]

barcodes = [f"{i:02}" for i in range(1, num_samples + 1)]
bin_size = 1000  # can be changed 

# group data by bin + mean depth in each bin 
def bin_coverage(df, bin_size):
    df = df.copy()
    startpos = df["pos"].min() 
    df["bin"] = ((df["pos"] - startpos) // bin_size) * bin_size + startpos 
    return df.groupby("bin")["depth"].mean().reset_index() 

coverage = []
bins = None
sample_labels = []

for barcode in barcodes:
    file_path = os.path.join(input, f"depth_barcode{barcode}.txt")
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        try:
            df = pd.read_csv(file_path, sep='\t', header=None, names=["chrom", "pos", "depth"])
            bf = bin_coverage(df, bin_size)
            if bins is None:
                bins = bf["bin"]
            coverage.append(binned_df["depth"].values)
            sample_labels.append(f"Sample {barcode}")
        except Exception as err:
            print(f"Couldn't read {file_path}: {err}")
    else:
        print(f"Non existing: {file_path}")

array = np.array(coverage)

plt.figure(figsize=(15, 6)) 
sns.heatmap(array, xticklabels=False, yticklabels=sample_labels, cmap="BuPu",cbar_kws={'label': 'Coverage depth'})
interval = 10  # can be changed --> which bins are shown on x axis
positions = np.arange(0, len(bins), interval)
labels = bin_positions[::interval].astype(str) 

plt.xticks(positions, labels, rotation=90)
plt.xlabel("Genomic position (1000 bp bins)") 
plt.title(f"Coverage heatmap: {args['region']}")
plt.tight_layout()
os.makedirs(output, exist_ok=True)
outputfile = os.path.join(output, f"coverage_heatmap_{region}.png")
plt.savefig(outputfile, dpi=300)
print(f"heatmap done ; {outputfile}")
