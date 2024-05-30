# HPC UGent ONT PreProcessing
This script is designed to perform various analyses on nanopore sequencing data, including basecalling, demultiplexing, alignment, quality control, and CNV (Copy Number Variation) analysis using tools such as Dorado, Minimap2, Samtools, PycoQC, Modkit, and WisecondorX.

## Prerequisites

Linux environment with Bash shell
SLURM job scheduler
Modules:
GCCcore/11.3.0
GCC/11.3.0
SAMtools/1.16.1-GCC-11.3.0
pod5-file-format/0.1.8-foss-2022a
minimap2/2.24-GCCcore-11.3.0

## Usage
Modify the script's parameters according to your specific dataset and analysis requirements.

./Main_Script_Test6-J.sh -p <InputDataPath> -t <InputDataType> -w <WORKDIR> -n <num_samples> -r <REF> -k <KIT_NAME> -c <CONFIG> -W <WISECONDORREF>

-p <InputDataPath>: Path to input data.
-t <InputDataType>: Type of input data (fast5 or pod5).
-w <WORKDIR>: Working directory where analysis results will be stored.
-n <num_samples>: Number of samples to analyze.
-r <REF>: Reference genome file.
-k <KIT_NAME>: Nanopore sequencing kit name (optional, default: SQK-NBD114-24).
-c <CONFIG>: Configuration file (optional, default: dna_r10.4.1_e8.2_400bps_sup@v4.2.0).
-W <WISECONDORREF>: Path to the WisecondorX reference file (optional, default: /kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/WisecondorX_ref/LQB.GRCh38.100kb.npz).

## Background Processes

The script utilizes SLURM to run background processes for basecalling and methylation analysis in parallel.

## Output

Basecalling results: Stored in ${WORKDIR}/basecalling.
Methylation analysis results: Stored in ${WORKDIR}/methylation.
Read length information: Stored in ${WORKDIR}/ReadLengths.
WisecondorX CNV analysis results: Stored in ${WORKDIR}/WisecondorX.

## Dependencies
Ensure the following dependencies are installed and available in your environment:

Dorado
Minimap2
Samtools
PycoQC
Modkit
WisecondorX



