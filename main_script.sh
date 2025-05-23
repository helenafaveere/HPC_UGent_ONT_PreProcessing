#!/bin/bash -l 

#SBATCH --time=72:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=24
#SBATCH --gpus-per-task=4
#SBATCH --mem-per-gpu=80G
#SBATCH --output=main_script_%j.log 
#SBATCH --error=main_script_%j.err 
#SBATCH --mail-user=helena.faveere@ugent.be
#SBATCH --mail-type=ALL
#SBATCH --job-name=main_script

###------------------------------------------------- module loading
module purge
module load GCCcore/12.3.0
module load GCC/12.3.0
module load SAMtools/1.18-GCC-12.3.0
module load pod5-file-format/0.3.10-foss-2023a
module load minimap2/2.26-GCCcore-12.3.0
module load R/4.3.3-gfbf-2023b

export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/dorado-0.8.2-linux-x64/bin:${PATH} # Add path to your Dorado installation
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/pycoqc_env/bin:${PATH} # Add path to your pycoQC installation
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/modkit_env/bin:${PATH} # Add path to your modkit installation
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/Rlibs/bin:${PATH} # Add path to your R environment
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/wisecondorx_env/bin:${PATH} # Add path to your WisecondorX installation
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/methylartist_env/bin:${PATH} # Add path to your Methylartist installation
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/nanocaller_env/bin:${PATH} # Add path to your Nanocaller installation
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/whatshap_env/bin:${PATH} # Add path to your Whatshap installation

###------------------------------------------------- flag definition
# Check for provided options
while getopts "p:t:w:n:r:k:c:W:" flag; do
    case "${flag}" in
        p) InputDataPath="${OPTARG}" ;;
        t) InputDataType="${OPTARG}" ;;
        w) WORKDIR="${OPTARG}" ;;
        n) num_samples="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        k) KIT_NAME="${OPTARG}" ;;
        c) CONFIG="${OPTARG}" ;;
        W) WISECONDORREF="${OPTARG}" ;; 
    esac
done

#default values
DEFAULT_KIT_NAME="SQK-NBD114-96"
DEFAULT_CONFIG="dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
DEFAULT_WISECONDORREF="/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/WisecondorX_ref/LQB.GRCh38.100kb.npz"

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
    pod5 convert fast5 -r -O "$InputDataPath" "$InputDataPath"/*.fast5 "$InputDataPath"/output_pod5s
    #Fast5_2_Pod5 "$InputDataPath"
    READSDIR="$InputDataPath"/"output_pod5s"
elif [ "$InputDataType" = "pod5" ]; then
    echo "Handling pod5 data type"
    READSDIR="$InputDataPath"
else
    echo "Unsupported data type: $InputDataType"
    exit 1
fi

###----------------------------------------------- exporting all variables to make them available inside srun
inside srun
export READSDIR=${READSDIR}
export CONFIG=${CONFIG}
export batch_size=${batch_size}
export chunk_size=${chunk_size}
export KIT_NAME=${KIT_NAME}
export num_samples=${num_samples}

export WORKDIR=${WORKDIR}
mkdir -p ${WORKDIR}
cd ${WORKDIR}
mkdir -p "${WORKDIR}/basecalling/demux/" 
mkdir -p "${WORKDIR}/methylation/demux/" 

###----------------------------------------------- settings for basecalling
echo ${CONFIG} 
echo ${KIT_NAME}
batch_size=0
chunk_size=10000
dorado download --model ${CONFIG}

###--------------------------------------------- simplex basecalling (background process)
srun --ntasks=1 --exclusive bash -c '
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
    dorado demux --output-dir "${WORKDIR}/basecalling/demux" --no-classify basecalling/simplex_all_barcodes.bam &&
    for ((i=1; i<=num_samples; i++)); do
        sample_num=$(printf "%02d" "$i")
        dorado summary -v "basecalling/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "basecalling/sequencing_summary_simplex_barcode${sample_num}.txt" 
    done
' &

###--------------------------------------------- methylation analysis (background process)
srun --ntasks=1 --exclusive bash -c '
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
    dorado demux --output-dir "${WORKDIR}/methylation/demux" --no-classify methylation/methylation_all_barcodes.bam &&
    for ((i=1; i<=num_samples; i++)); do
        sample_num=$(printf "%02d" "$i")
        dorado summary -v "methylation/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "methylation/sequencing_summary_methylation_barcode${sample_num}.txt"
    done
' &
wait

###------------------------------------------------ settings for basecalling alignment
export REF=${REF}

###--------------------------------------------- alignment and QC for basecalling
for ((i=1; i<=num_samples; i++)); do
    sample_num=$(printf "%02d" "$i")
    
    samtools bam2fq "basecalling/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "basecalling/${KIT_NAME}_barcode${sample_num}.fastq"

    dorado aligner "${REF}" "basecalling/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "basecalling/simplex_mapped_barcode${sample_num}.bam"

    samtools sort "basecalling/simplex_mapped_barcode${sample_num}.bam" -o "basecalling/sorted_simplex_mapped_barcode${sample_num}.bam"
    
    samtools index "basecalling/sorted_simplex_mapped_barcode${sample_num}.bam"

    pycoQC -f "basecalling/sequencing_summary_simplex_barcode${sample_num}.txt" -a "basecalling/sorted_simplex_mapped_barcode${sample_num}.bam" -o "basecalling/QC_simplex_barcode${sample_num}.html"
done

###------------------------------------------------ settings for methylation alignment
export REF=${REF}

###--------------------------------------------- alignment for methylation analysis
for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    
    samtools bam2fq "methylation/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "methylation/${KIT_NAME}_barcode${sample_num}.fastq"

    dorado aligner "${REF}" "methylation/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "methylation/mapped_methylation_barcode${sample_num}.bam"
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
mkdir -p "${WORKDIR}/ReadLengths" 
for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    
    samtools view "basecalling/${KIT_NAME}_barcode${sample_num}.bam" | awk '{print length($10)}' > "ReadLenghts/all_read_lengths_barcode${sample_num}.txt"
    samtools view "basecalling/simplex_mapped_barcode${sample_num}.bam" | awk '{print length($10)}' > "ReadLenghts/all_mapped_read_lengths_barcode${sample_num}.txt"
done


###----------------------------------------------- running WisecondorX to study CNVs (generates plots and BED files)
 # conda activate $CONDA_ENV_NAME #activate wisecondorx environment # now at the top

export WISECONDORREF=${WISECONDORREF}
mkdir -p "${WORKDIR}/WisecondorX" 

for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    output_folder="barcode${i}"
    mkdir -p "${WORKDIR}/WisecondorX/${output_folder}" 
    
    WisecondorX convert "basecalling/simplex_mapped_barcode${sample_num}.bam" "WisecondorX/simplex_mapped_barcode${sample_num}.npz"
    
    WisecondorX predict "WisecondorX/simplex_mapped_barcode${sample_num}.npz" ${WISECONDORREF} barcode${sample_num} --plot --bed
    
    mv barcode${sample_num}* "${WORKDIR}/WisecondorX/${output_folder}/"

done
