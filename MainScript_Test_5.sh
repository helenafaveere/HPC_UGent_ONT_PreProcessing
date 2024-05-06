#!/bin/bash -l

#SBATCH --time=72:00:00
#SBATCH --mem=250gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16   

#SBATCH --gres=gpu:4
#SBATCH --mail-user=justine.rayp@ugent.be
#SBATCH --mail-type=ALL
#SBATCH --job-name=ONT-Test5-J

###------------------------------------------------- module loading
module purge
module unload GCC/13.2.0
module unload GCCcore/13.2.0
module unload zlib/1.2.13-GCCcore-13.2.0
module unload binutils/2.40-GCCcore-13.2.0
module load GCCcore/11.2.0
module load GCC/11.2.0
module load SAMtools/1.14-GCC-11.2.0
module load minimap2/2.22-GCCcore-11.2.0
export PATH=/scratch/gent/vo/001/gvo00115/vsc45900/dorado-0.4.0-linux-x64/bin:${PATH} # Add path to dorado installation
export PATH=/scratch/gent/vo/001/gvo00115/vsc45900/miniforge3/envs/pycoqc_252/bin:${PATH} # Set this to the right path to your installation of pycoqc
export PATH=/scratch/gent/vo/001/gvo00115/vsc45900/miniforge3/envs/modkit/bin:${PATH} # Add the path to your modkit installation
conda activate wisecodorXv1.2.5 #activate your wisecondorx environment

###------------------------------------------------- function definition
Fast5_2_Pod5() {
   # Convert Pod5 to Fast5
   pod5 convert fast5 "${1}"/*.fast5 output_pod5s --one-to-one ./
}

###------------------------------------------------- flag definition
DEFAULT_KIT_NAME="SQK-NBD114-24"
DEFAULT_CONFIG="dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
DEFAULT_WISECONDORREF="/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/WisecondorX_ref/LQB.GRCh38.100kb.npz"

# Check for provided options
while getopts "p:t:w:n:r:k:c:W:" flag; do
    case "${flag}" in
        p) InputDataPath="${OPTARG}" ;;
        t) InputDataType="${OPTARG}" ;;
        w) WORKDIR="${OPTARG}" ;;
        n) num_samples="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        k) KIT_NAME="${OPTARG:-"SQK-NBD114-24"}" ;;
        c) CONFIG="${OPTARG:-"dna_r10.4.1_e8.2_400bps_sup@v4.2.0"}" ;;
        W) WISECONDORREF="${OPTARG:-"/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/WisecondorX_ref/LQB.GRCh38.100kb.npz"}" ;; 
    esac
done

# Set default values if not provided
if [[ -z "$KIT_NAME" ]]; then
    KIT_NAME="$DEFAULT_KIT_NAME"
fi

if [[ -z "$CONFIG" ]]; then
    CONFIG="$DEFAULT_CONFIG"
fi
if [[ -z "$WISECONDORREF" ]]; then
    WISECONDORREF="$DEFAULT_WISECONDORREF"
fi

echo "KIT_NAME: $KIT_NAME"
echo "CONFIG: $CONFIG"
echo "WISECONDORREF: $WISECONDORREF"


###------------------------------------------------- converting to the right datatype
if [ "$InputDataType" = "fast5" ]; then
    Fast5_2_Pod5 "$InputDataPath"
    READSDIR="/outputpod5s"
elif [ "$InputDataType" = "pod5" ]; then
    echo "Handling pod5 data type"
    READSDIR="$InputDataPath"
else
    echo "Unsupported data type: $InputDataType"
    exit 1
fi


###----------------------------------------------- exporting all variables to make them available inside srun

export READSDIR=${READSDIR}
export CONFIG=${CONFIG}
export batch_size=${batch_size}
export chunk_size=${chunk_size}
export KIT_NAME=${KIT_NAME}
export num_samples=${num_samples}

export WORKDIR=${WORKDIR}
mkdir -p ${WORKDIR}
cd ${WORKDIR}
mkdir -p "${WORKDIR}/basecalling" 
mkdir -p "${WORKDIR}/methylation" 

###----------------------------------------------- settings for basecalling
echo ${CONFIG} 
echo ${KIT_NAME}
batch_size=0
chunk_size=10000
dorado download --model ${CONFIG}

###--------------------------------------------- simplex basecalling (background process)

srun --cpu-bind=cores --cpus-per-task=16 bash -c '
    echo "start basecalling" &&
    echo ${CONFIG} &&
    echo ${KIT_NAME} &&
    echo ${WORKDIR}/${CONFIG} &&
    echo ${READSDIR} &&
    dorado basecaller \
        --device "cuda:0,1,2,3" \
        --batchsize ${batch_size} \
        --chunksize ${chunk_size} \
        --recursive \
        --kit-name ${KIT_NAME} \
        --verbose \
        ${WORKDIR}/${CONFIG} \
        ${READSDIR} > basecalling/simplex_all_barcodes.bam &&
    dorado demux --kit-name "${KIT_NAME}" --output-dir "${WORKDIR}/basecalling" basecalling/simplex_all_barcodes.bam &&
    for ((i=1; i<=num_samples; i++)); do
        sample_num=$(printf "%02d" "$i")
        dorado summary -v "basecalling/${KIT_NAME}_barcode${sample_num}.bam" > "basecalling/sequencing_summary_simplex_barcode${sample_num}.txt" 
    done
' &

###--------------------------------------------- methylation analysis (background process)

srun --cpu-bind=cores --cpus-per-task=16 bash -c '
    dorado basecaller \
        --device "cuda:0,1,2,3" \
        --batchsize "${batch_size}" \
        --chunksize "${chunk_size}" \
        --modified-bases 5mCG_5hmCG \
        --recursive \
        --kit-name "${KIT_NAME}" \
        --verbose \
        "${WORKDIR}/${CONFIG}" \
        "${READSDIR}" > methylation/methylation_all_barcodes.bam &&
    dorado demux --kit-name "${KIT_NAME}" --output-dir "${WORKDIR}/methylation" methylation/methylation_all_barcodes.bam &&
    for ((i=1; i<=num_samples; i++)); do
        sample_num=$(printf "%02d" "$i")
        dorado summary -v "methylation/${KIT_NAME}_barcode${sample_num}.bam" > "methylation/sequencing_summary_methylation_barcode${sample_num}.txt"
    done
' &

wait

###------------------------------------------------ settings for basecalling alignment

export REF=${REF}

###--------------------------------------------- alignment and QC for basecalling

for ((i=1; i<=num_samples; i++)); do
    sample_num=$(printf "%02d" "$i")
    
    samtools bam2fq "basecalling/${KIT_NAME}_barcode${sample_num}.bam" > "basecalling/${KIT_NAME}_barcode${sample_num}.fastq"

    minimap2 -ax map-ont -t 30 "${REF}" "basecalling/${KIT_NAME}_barcode${sample_num}.fastq" | samtools sort -o "basecalling/simplex_mapped_barcode${sample_num}.bam"
    
    samtools index "basecalling/simplex_mapped_barcode${sample_num}.bam"

    pycoQC -f "basecalling/sequencing_summary_simplex_barcode${sample_num}.txt" -a "basecalling/simplex_mapped_barcode${sample_num}.bam" -o "basecalling/QC_simplex_barcode${sample_num}.html"
done

###------------------------------------------------ settings for methylation alignment

###--------------------------------------------- alignment for methylation analysis

for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    
    samtools bam2fq "methylation/${KIT_NAME}_barcode${sample_num}.bam" > "methylation/${KIT_NAME}_barcode${sample_num}.fastq"

    dorado aligner "${REF}" "methylation/${KIT_NAME}_barcode${sample_num}.bam" > "methylation/mapped_methylation_barcode${sample_num}.bam"
done

###--------------------------------------------- methylation QC

for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    
    samtools sort "methylation/mapped_methylation_barcode${sample_num}.bam" -o "methylation/sorted_methylation_barcode${sample_num}.bam"
    samtools index "methylation/sorted_methylation_barcode${sample_num}.bam"
    
    modkit pileup "methylation/sorted_methylation_barcode${sample_num}.bam" "methylation/methylation_barcode${sample_num}.bed" \
        --ref "${REF}" \
        --preset traditional
done

###----------------------------------------------- get the readlenghts for further processing
for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    
    samtools view SQK-NBD114-24_barcode${sample_num}.bam | awk '{print length($10)}' > all_read_lengths_barcode${sample_num}.txt
done


###----------------------------------------------- running WisecondorX to study CNVs (generates plots and BED files)

export WISECONDORREF=${WISECONDORREF}


for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    
    WisecondorX convert simplex_mapped_barcode${sample_num}.bam simplex_mapped_barcode${sample_num}.npz
    
    WisecondorX predict simplex_mapped_barcode${sample_num}.npz ${WISECONDORREF} barcode${sample_num} --plot --bed
    
done
