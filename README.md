# UGent_cfDNA_ONT_pipeline

This pipeline is adapted from *HPC UGent ONT PreProcessing by JustineRayp* and is developed to perform various analyses on ONT sequencing data, like basecalling, demultiplexing, alignment, quality control, copy number variation analysis, and methylation analysis. A range of tools are employed, including dorado, modkit, samtools, pycoQC, QDNAseq, wisecondorX, PREFACE, nanocaller, whatshap, and methylartist.

![Image](https://github.com/user-attachments/assets/b1a8e15d-a2b5-4743-8939-c824884de39d)

## Prerequisites

- Linux environment with Bash shell
- SLURM job scheduler
- Modules:
  - GCCcore/13.3.0
  - GCC/13.3.0 
  - SAMtools/1.21-GCC-13.3.0 
  - pod5-file-format/0.3.10-foss-2023a
  - R/4.4.2-gfbf-2024a

## Usage
Modify the script's parameters according to your specific dataset and analysis requirements.
```
./main_script.sh -p -t -w -n <num_samples> -r -k <KIT_NAME> -c -W -g -h
```

- -p <InputDataPath>: Path to input data.
- -t <InputDataType>: Type of input data (FAST5 or POD5).
- -w <WORKDIR>: Working directory where analysis results will be stored.
- -n <num_samples>: Number of samples to analyze.
- -r <REF>: Reference genome file.
- -k <KIT_NAME>: Nanopore sequencing kit name (optional, default: SQK-NBD114-24).
- -c <CONFIG>: Configuration file (optional, default: dna_r10.4.1_e8.2_400bps_sup@v4.2.0).
- -W <WISECONDORREF>: Path to the WisecondorX reference file (optional, default: /kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/WisecondorX_ref/LQB.GRCh38.100kb.npz).
- -g <REGION>: Gene region (optional, default: chrX:73820672-73838873)
- -h <GENE>: Gene name corresponding to the region (optional, default: XIST)

## Background Processes

The script utilizes SLURM to run background processes for basecalling and methylation analysis in parallel.

## Output

- Basecalling results: Stored in `${WORKDIR}/basecalling`.
- Methylation analysis results: Stored in `${WORKDIR}/methylation`.
- Read length information: Stored in `${WORKDIR}/ReadLengths`.
- Histograms of the read lengths: Stored in `${WORKDIR}/ReadLengths/Histogram`.
- WisecondorX CNV analysis results: Stored in `${WORKDIR}/WisecondorX`.
- QDNAseq CNV analysis results: Stored in `${WORKDIR}/QDNAseq`.
- PREFACE results: Stored in `${WORKDIR}/PREFACE`.
- Heatmap of the coverage within gene region: Stored in `${WORKDIR}/Heatmaps`.
- Results nanocaller: Stored in `${WORKDIR}/NanoCaller`.
- Haplotagged BAM files: Stored in `${WORKDIR}/WhatsHap`.
- Output methylartist locus: Stored in `${WORKDIR}/Methylartist`.

## Dependencies
Ensure the following dependencies are installed and available in your (conda) environment:

- Dorado
- Samtools
- PycoQC
- Modkit
- WisecondorX
- Python
- R
- Methylartist
- Nanocaller
- Whatshap



