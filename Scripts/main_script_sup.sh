#!/bin/bash -l 

#SBATCH --time=72:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=24
#SBATCH --gpus-per-task=4
#SBATCH --mem-per-gpu=80G
#SBATCH --output=TestRun_sup_%j.log # Output log file (%j will be replaced by the job ID)
#SBATCH --error=TestRun_sup_%j.err  # Error log file (%j will be replaced by the job ID)
#SBATCH --mail-user=#fillin@ugent.be
#SBATCH --mail-type=ALL
#SBATCH --job-name=SUP

###------------------------------------------------- module loading
module purge
module load GCCcore/13.3.0
module load GCC/13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load pod5-file-format/0.3.10-foss-2023a
module load minimap2/2.26-GCCcore-13.3.0
module load R/4.4.2-gfbf-2024a

export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/dorado-0.8.2-linux-x64/bin:${PATH} # Add path to your Dorado installation

###------------------------------------------------- flag definition
# Check for provided options
while getopts "p:t:w:n:r:k:c:W:g:h:" flag; do
    case "${flag}" in
        p) InputDataPath="${OPTARG}" ;;
        t) InputDataType="${OPTARG}" ;;
        w) WORKDIR="${OPTARG}" ;;
        n) num_samples="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        k) KIT_NAME="${OPTARG}" ;;
        c) CONFIG="${OPTARG}" ;;
        W) WISECONDORREF="${OPTARG}" ;; 
        g) REGION="${OPTARG}" ;; 
        h) GENE="${OPTARG}" ;; 
    esac
done

#default values
DEFAULT_KIT_NAME="SQK-NBD114-24"
DEFAULT_CONFIG="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
DEFAULT_WISECONDORREF="/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/WisecondorX_ref/LQB.GRCh38.100kb.npz" # Change to the location of your WisecondorX reference
DEFAULT_REGION="chrX:73820672-73838873" 
DEFAULT_GENE="XIST"

if [[ -z "$KIT_NAME" ]]; then
    KIT_NAME="$DEFAULT_KIT_NAME"
fi

if [[ -z "$CONFIG" ]]; then
    CONFIG="$DEFAULT_CONFIG"
fi
if [[ -z "$WISECONDORREF" ]]; then
    WISECONDORREF="$DEFAULT_WISECONDORREF"
fi
if [[ -z "$REGION" ]]; then
    REGION="$DEFAULT_REGION"
fi
if [[ -z "$GENE" ]]; then
    GENE="$DEFAULT_GENE"
fi
echo "KIT_NAME: $KIT_NAME"
echo "CONFIG: $CONFIG"
echo "WISECONDORREF: $WISECONDORREF"
echo "REGION: $REGION"
echo "GENE: $GENE"

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
echo "demux + basecall + summary for both types done"

###------------------------------------------------ settings for basecalling alignment
export REF=${REF}
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/pycoqc_env/bin:${PATH} # Add path to your pycoQC installation

###--------------------------------------------- alignment and QC for basecalling
for ((i=1; i<=num_samples; i++)); do
    sample_num=$(printf "%02d" "$i")
    
    samtools bam2fq "basecalling/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "basecalling/${KIT_NAME}_barcode${sample_num}.fastq"
    
    dorado aligner "${REF}" "basecalling/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "basecalling/simplex_mapped_barcode${sample_num}.bam"

    samtools sort "basecalling/simplex_mapped_barcode${sample_num}.bam" -o "basecalling/sorted_simplex_mapped_barcode${sample_num}.bam"
    
    samtools index "basecalling/sorted_simplex_mapped_barcode${sample_num}.bam"

    pycoQC -f "basecalling/sequencing_summary_simplex_barcode${sample_num}.txt" -a "basecalling/sorted_simplex_mapped_barcode${sample_num}.bam" -o "basecalling/QC_simplex_barcode${sample_num}.html"
done

echo "QC basecalling ok"
###------------------------------------------------ settings for methylation alignment
export REF=${REF}

###--------------------------------------------- alignment for methylation analysis
for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    
    samtools bam2fq "methylation/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "methylation/${KIT_NAME}_barcode${sample_num}.fastq"

    dorado aligner "${REF}" "methylation/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" > "methylation/mapped_methylation_barcode${sample_num}.bam"
done

echo "methylation alignment ok"

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

echo "methylation QC ok"

###----------------------------------------------- get the readlenghts for further processing
mkdir -p "${WORKDIR}/ReadLengths" 
for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    
    samtools view "basecalling/demux/5b6469e0391acf348cd89728e70975aabd01996f_${KIT_NAME}_barcode${sample_num}.bam" | awk '{print length($10)}' > "ReadLengths/all_read_lengths_barcode${sample_num}.txt"
    samtools view "basecalling/sorted_simplex_mapped_barcode${sample_num}.bam" | awk '{print length($10)}' > "ReadLengths/all_mapped_read_lengths_barcode${sample_num}.txt"
done

echo "readlengths ok"

###----------------------------------------------- make histograms of the readlengths 
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/python_env/bin:${PATH} # Add path to your python environment

mkdir -p "${WORKDIR}/ReadLengths/Histogram"

for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    python /kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Scripts/generate_read_lengths_histograms.py "$sample_num" "${WORKDIR}"
done

echo "histograms ok"

###----------------------------------------------- running WisecondorX to study CNVs (generates plots and BED files)
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/wisecondorx_env/bin:${PATH} # Add path to your WisecondorX installation

export WISECONDORREF=${WISECONDORREF}
mkdir -p "${WORKDIR}/WisecondorX" 

for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    output_folder="barcode${i}"
    mkdir -p "${WORKDIR}/WisecondorX/${output_folder}" 
    
    WisecondorX convert "basecalling/sorted_simplex_mapped_barcode${sample_num}.bam" "WisecondorX/simplex_mapped_barcode${sample_num}.npz"
    
    WisecondorX predict "WisecondorX/simplex_mapped_barcode${sample_num}.npz" ${WISECONDORREF} barcode${sample_num} --plot --bed
    
    mv barcode${sample_num}* "${WORKDIR}/WisecondorX/${output_folder}/"

done

echo "wisecondorx ok"
###----------------------------------------------- running QDNAseq to study CNVs
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/Rlibs/bin:${PATH} # Add path to your R environment

mkdir -p "${WORKDIR}/QDNAseq" 

for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    Rscript /kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Scripts/GRCh38_VisualisationScript.R "$WORKDIR" "$sample_num"
done

echo "qdnaseq ok"

###----------------------------------------------- running PREFACE to study FF
mkdir -p "${WORKDIR}/PREFACE/bedfiles"

for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "barcode%02d" "$i")
    if [[ -f "$infile" ]]; then
        awk 'BEGIN {OFS="\t"} {gsub(/nan/, "NA", $5); print $1, $2, $3, $5}' "${WORKDIR}/WisecondorX/${barcode}/${sample_num}_bins.bed" > "${WORKDIR}/PREFACE/bedfiles/${sample_num}_fixed.bed"
    else
        echo "Warning: $infile not found." >&2
    fi
done

for bedfile in "${WORKDIR}/PREFACE/bedfiles/"*_fixed.bed; do
    sample_name=$(basename "$bedfile" _fixed.bed)
    out_prefix="${WORKDIR}/PREFACE/${sample_name}"
    echo "Running PREFACE on $sample_name"
    Rscript "/user/gent/484/vsc48405/gONT/Scripts/PREFACE.R" predict --infile "$bedfile" --model "/user/gent/484/vsc48405/gONT/Scripts/PREFACE.RData" --json "$out_prefix"
done

echo "preface ok"

###----------------------------------------------- running methylation heatmap
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/python_env/bin:${PATH} # Add path to your python environment

mkdir -p "${WORKDIR}/Heatmaps/coverage_outputs_samtools"

for ((i=1; i<=num_samples; i++))
do
    sample_num=$(printf "%02d" "$i")
    bam_file="methylation/sorted_methylation_barcode${sample_num}.bam"
    output="Heatmaps/coverage_outputs_samtools/depth_barcode${sample_num}.txt"

    if [ -f "${bam_file}" ]; then
        echo "Calculating depth for ${bam_file} in region ${REGION}..."
        samtools depth -r "${REGION}" -a "${bam_file}" > "${output}"
        echo "Saved depth output to ${output}"
    else
        echo "Warning: ${bam_file} not found, skipping..."
    fi
done

python "/user/gent/484/vsc48405/gONT/Scripts/generate_heatmap.py" \
    -i "${WORKDIR}/Heatmaps/coverage_outputs_samtools" \
    -r "${REGION}" \
    -o "${WORKDIR}/Heatmaps" \
    -n "${num_samples}"
    
echo "heatmap ok"

###----------------------------------------------- running NanoCaller
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/nanocaller_env/bin:$PATH
export REF=${REF}

mkdir -p "${WORKDIR}/NanoCaller"

for ((i=1; i<=num_samples; i++)); do
    sample_num=$(printf "%02d" "$i")

    NanoCaller --bam "${WORKDIR}/methylation/sorted_methylation_barcode${sample_num}.bam" \
               --ref "$REF" \
               --preset ont \
               --phase \
               --enable_whatshap \
               --sequencing short_ont \
               --mincov 1 \
               --regions "$REGION" \
               --output "${WORKDIR}/NanoCaller/barcode${sample_num}.vcf" \
               --prefix "barcode${sample_num}"
done

echo "nanocaller ok"

###----------------------------------------------- running whatshap haplotag
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/whatshap_env/bin:$PATH
export REF=${REF}

mkdir -p "${WORKDIR}/WhatsHap"

for ((i=1; i<=num_samples; i++)); do
    sample_num=$(printf "%02d" "$i")
    
    whatshap haplotag \
        -o "${WORKDIR}/WhatsHap/haplotagged_barcode${sample_num}.bam" \
        --reference "$REF" \
        "${WORKDIR}/NanoCaller/barcode${sample_num}.vcf/barcode${sample_num}.snps.phased.vcf.gz" \
        "${WORKDIR}/methylation/sorted_methylation_barcode${sample_num}.bam" \
        --ignore-read-groups \
        --skip-missing-contigs 
        
    samtools index "${WORKDIR}/WhatsHap/haplotagged_barcode${sample_num}.bam"
done

echo "whatshap haplotag ok"

###----------------------------------------------- running methylartist 
export PATH=/kyukon/data/gent/shared/001/gvo00115/ONT_cfDNA/Tools/miniforge3/envs/methylartist_env/bin:$PATH

mkdir -p "${WORKDIR}/Methylartist"

for ((i=1; i<=num_samples; i++)); do
    sample_num=$(printf "%02d" "$i")
    
    methylartist locus \
            -b "${WORKDIR}/WhatsHap/haplotagged_barcode${sample_num}.bam" \
            -i "$REGION" \
            --gtf "/user/gent/484/vsc48405/gONT/Tools/ncbi_dataset/hg38_annotation/Genes.gencode/sorted.gtf.gz" \
            --genes "$GENE" \
            --labelgenes \
            -p 1,6,1,3,4 \
            --ref "$REF" \
            --motif CG \
            --phased \
            -o "${WORKDIR}/Methylartist/methylartist_${GENE}_${sample_num}"
done

echo "methylartist ok"
