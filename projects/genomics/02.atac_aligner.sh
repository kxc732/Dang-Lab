#!/bin/bash

# FILE INFORMATION
#---------------------------
# Author: Kasonde Chewe
# Pipeline: ATAC-seq 0 2 24, knockdown, FL, DLCS1, DLCS2
# File: 02.atac_aligner.sh 
# Date: 05/01/2024
# Last Modification: 05/12/2024
# Source Code Location: ~/final-project-dir/ATAC-SEQ/scripts/ultimate/bin
# Description: Process RAW ATAC-SEQ from FASTQ files to generate BAM and counts files
# Additional files generated 


# SLURM DIRECTIVES NODE06 working fine 
#---------------------------
#SBATCH --job-name=chromatin_accessibility
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=40  # Since you have 40 tasks per node
#SBATCH --mem=128G   # Memory requirement


# Quality Control 
# --------------------------------
# 1. fastqc verified that adapter seuqences where from NexteraPE-PE.fa 
# 2. set minimum quality score phred 20 and minimum length 50 bp 
# 3. picard mark duplicates 
# 4. bedtools intersect and filter reads 
# 5. samtools flags etc 
# 6. macs2 vs macs3 + params 

 
# [TO-DO-1]: add IDR quality control 



# PARSE COMMAND LINE ARGUMENTS
#---------------------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --list)
            IFS=',' read -r -a file_list <<< "${2//\'/}"
            shift
            ;;
        --fastqc_out_prefix) fastqc_out_prefix="${2//\'/}"; shift ;;
        --trimmed_out_prefix) trimmed_out_prefix="${2//\'/}"; shift ;;
        --bam_out_prefix) bam_out_prefix="${2//\'/}"; shift ;;
        *)
            echo "Unknown parameter passed: $1"
            exit 1
            ;;
    esac
    shift
done


# Make logging directory in script folder
# ----------------------------------------
mkdir -p ./logs

# ACTIVATE CONDA ENVIRONMENT
#----------------------------
eval "$(conda shell.bash hook)"
conda activate gro-seq

# Initialize counters
#---------------------------
total_samples=0
successful_alignments=0
missing_pairs=0

# VARIABLES
#---------------------------
DIR=~/final-project-dir/ATAC-SEQ/ATACSEQ
ADAPTER=~/final-project-dir/ATAC-SEQ/scripts/ultimate/utilities/adapters/NexteraPE-PE.fa
BLACKLIST=~/final-project-dir/ATAC-SEQ/scripts/ultimate/utilities/blacklist_regions/RegionsRemove.bed
FASTQC_DIR=${fastqc_out_prefix%\"}
TRIMMED_DIR=${trimmed_out_prefix%\"}
POST_FASTQC_DIR=$DIR/fastqc/post
GENOME_DIR=~/genomes/hg38
INDEX_DIR=~/genomes/hg38/index/bwa/hg38.fa
GENOME_FILE=$GENOME_DIR/hg38.fa
GTF_FILE=$GENOME_DIR/annotations/hg38.ncbiRefSeq.gtf
BAMS_DIR=${bam_out_prefix%\"}
THREADS=40



mkdir -p $FASTQC_DIR $TRIMMED_DIR $BAMS_DIR $POST_FASTQC_DIR $BAMS_DIR/samtools_stats

# LOGGING
#----------------------------
LOG_DIR=$DIR/logs
mkdir -p $LOG_DIR $LOG_DIR/report
LOG=$LOG_DIR/log_atac_seq_aligner_$(date +%Y-%m-%d%-H:%M:%S).txt
touch $LOG
ulimit -n 4096

# PROCESSING
#----------------------------
echo -e "\n Starting ATAC-seq processing pipeline. [UPGRADE-3] \n ======================= \n" | tee -a $LOG

# Sample processing loop just for testing file pairing
for ((i=0; i<${#file_list[@]}; i+=2)); do
    FASTQ_R1=${file_list[i]}
    FASTQ_R2=${file_list[i+1]}

    if [[ -z "$FASTQ_R1" || -z "$FASTQ_R2" ]]; then
        echo -e "[ERROR] -- $(date +%H:%M:%S) --: Missing pair for sample ${SAMPLE}. R1 or R2 file not found." | tee -a $LOG
        ((missing_pairs++))
        continue
    else
        R1=$(basename $FASTQ_R1 "_L00[0-9]_R1_001.fastq.gz")
        R2=$(basename $FASTQ_R2 "_L00[0-9]_R2_001.fastq.gz")
        echo "[LOG] -- $(date +%H:%M:%S) -- Processing: $R1 and $R2" | tee -a $LOG
        ((total_samples++))
    fi

    # Running FastQC before trimming
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Running FastQC {skipped}..." | tee -a $LOG
    fastqc -t 2 -o $FASTQC_DIR --noextract $FASTQ_R1 $FASTQ_R2 2>> $LOG

    # Running Trimmomatic for read trimming
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Running fastp {skipped}..." | tee -a $LOG
    fastp -i $FASTQ_R1 -I $FASTQ_R2 -o $TRIMMED_DIR/${R1}_R1_paired.fastq.gz \
      -O $TRIMMED_DIR/${R2}_R2_paired.fastq.gz \
      --cut_front --cut_tail --cut_window_size 4 --trim_poly_g \
      --cut_mean_quality 20 --length_required 45 \
      --thread $THREADS 2>> $LOG


    # Running FastQC after trimming
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Running post-trimming FastQC {skipped}..." | tee -a $LOG
    fastqc -t 2 -o $POST_FASTQC_DIR --noextract \
           $TRIMMED_DIR/${R1}_R1_paired.fastq.gz \
           $TRIMMED_DIR/${R2}_R2_paired.fastq.gz 2>> $LOG


    # Running BWA for alignment
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Running BWA alignment..." | tee -a $LOG
    bwa mem -t $THREADS $INDEX_DIR $TRIMMED_DIR/${R1}_R1_paired.fastq.gz $TRIMMED_DIR/${R2}_R2_paired.fastq.gz \
        | samtools view -Sb - \
        | samtools sort -o $BAMS_DIR/${R1}_sorted.bam -
    samtools index $BAMS_DIR/${R1}_sorted.bam 2>> $LOG
    ((successful_alignments++))

    # Remove Mitochondrial reads
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Removing mitochondrial reads..." | tee -a $LOG
    samtools view -@ $THREADS -h $BAMS_DIR/${R1}_sorted.bam | grep -v chrM | \
        samtools sort -@ $THREADS -O bam -o $BAMS_DIR/${R1}_sorted.rmChrM.bam
    samtools index $BAMS_DIR/${R1}_sorted.rmChrM.bam 2>> $LOG

    # Mark duplicates with Picard on mitochondrial-removed BAM
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Marking duplicates on mitochondrial-removed BAM..." | tee -a $LOG
    picard MarkDuplicates I=$BAMS_DIR/${R1}_sorted.rmChrM.bam \
                          O=$BAMS_DIR/${R1}_sorted.rmChrM.dups.bam \
                          M=$BAMS_DIR/${R1}_sorted.rmChrM.metrics.txt \
                          REMOVE_DUPLICATES=true 2>> $LOG
    samtools index $BAMS_DIR/${R1}_sorted.rmChrM.dups.bam 2>> $LOG

    # Remove reads from blacklisted regions
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Removing reads from blacklisted regions..." | tee -a $LOG
    bedtools intersect -a $BAMS_DIR/${R1}_sorted.rmChrM.dups.bam -b $BLACKLIST -v -ubam > $BAMS_DIR/${R1}_final.bam
    samtools index $BAMS_DIR/${R1}_final.bam 2>> $LOG

    # Remove multi-mapped reads with MAPQ < 30
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Removing multi-mapped reads (MAPQ < 30)..." | tee -a $LOG
    samtools view -@ $THREADS -h -q 30 $BAMS_DIR/${R1}_final.bam > $BAMS_DIR/${R1}_final.rmMulti.bam
    samtools sort -@ $THREADS -o $BAMS_DIR/${R1}_final.rmMulti.sorted.bam $BAMS_DIR/${R1}_final.rmMulti.bam
    samtools index $BAMS_DIR/${R1}_final.rmMulti.sorted.bam 2>> $LOG
    #
    # Shift alignments using alignmentSieve for ATAC-seq data
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Shifting ATAC-seq alignments..." | tee -a $LOG
    alignmentSieve --numberOfProcessors $THREADS --ATACshift --bam $BAMS_DIR/${R1}_final.rmMulti.sorted.bam -o $BAMS_DIR/${R1}_shifted.tmp.bam

    # Sort the shifted BAM file and index it
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Sorting shifted BAM file..." | tee -a $LOG
    samtools sort -@ $THREADS -O bam -o $BAMS_DIR/${R1}_shifted.bam $BAMS_DIR/${R1}_shifted.tmp.bam
    samtools index -@ $THREADS $BAMS_DIR/${R1}_shifted.bam 2>> $LOG

    # Clean up temporary files
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Cleaning up temporary files..." | tee -a $LOG
    rm $BAMS_DIR/${R1}_shifted.tmp.bam

    # Flagstats for BAM statistics
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Generating BAM file statistics..." | tee -a $LOG
    samtools flagstat $BAMS_DIR/${R1}_sorted.bam > $BAMS_DIR/samtools_stats/${R1}_sorted.bam.stats
    samtools flagstat $BAMS_DIR/${R1}_sorted.rmChrM.dups.bam > $BAMS_DIR/samtools_stats/${R1}_sorted.rmChrM.dups.bam.stats
    samtools flagstat $BAMS_DIR/${R1}_final.rmMulti.sorted.bam > $BAMS_DIR/samtools_stats/${R1}_final.rmMulti.sorted.bam.stats
    samtools flagstat $BAMS_DIR/${R1}_shifted.bam > $BAMS_DIR/samtools_stats/${R1}_shifted.bam.stats

    # Running Count Matrix Generation
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Running HOMER makeTagDirectory..." | tee -a $LOG
    makeTagDirectory $BAMS_DIR/${R1}_tags $BAMS_DIR/${R1}_shifted.bam -genome $GENOME_DIR -checkGC 2>> $LOG

    echo "[LOG] -- $(date +%H:%M:%S) -- Generating counts matrix with HOMER..." | tee -a $LOG
    analyzeRepeats.pl hg38 -count exons -d $BAMS_DIR/${R1}_tags/ > $BAMS_DIR/${R1}_homerCounts.txt 2>> $LOG

    echo "[LOG] -- $(date +%H:%M:%S) -- Running featureCounts for peak regions..." | tee -a $LOG
    featureCounts -p -T $THREADS -a $GTF_FILE \
                -o $BAMS_DIR/${R1}_featureCounts.txt $BAMS_DIR/${R1}_shifted.bam 2>> $LOG

    # Running Tracks (Bigwig and Bedgraph Generation)
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Running HOMER makeUCSCfile to generate .bedgraph..." | tee -a $LOG
    makeUCSCfile $BAMS_DIR/${R1}_tags -o $BAMS_DIR/${R1}_track.bedGraph 2>> $LOG

    echo "[LOG] -- $(date +%H:%M:%S) -- Running deeptools bamCoverage to generate .BigWig..." | tee -a $LOG
    bamCoverage -b $BAMS_DIR/${R1}_final.rmMulti.sorted.bam  -o $BAMS_DIR/${R1}.bw --normalizeUsing RPKM -p $THREADS 2>> $LOG

    # Running PeakCaller (Bigwig and Bedgraph Generation)
    # --------------------------------------------------------------------------
    echo "[LOG] -- $(date +%H:%M:%S) -- Running peak caller MACS2 to generate .NarrowPeaks..." | tee -a $LOG
    macs2 callpeak -f BAMPE -g hs --keep-dup all --cutoff-analysis -n ${R1} -q 0.05 \
                   -t $BAMS_DIR/${R1}_shifted.bam --outdir $BAMS_DIR/${R1} 2>> $LOG
                   
    # Running IDR for peaks files 

done


# Final report
# --------------------------------------------------------------------------
echo -e "[REPORT]:Total missing pairs: $missing_pairs " | tee -a $LOG
echo -e "[REPORT]:Total samples processes: $total_samples " | tee -a $LOG
echo -e "[REPORT]:Total samples processes: $successful_alignments" | tee -a $LOG
echo "[FINISHED] -- $(date +%H:%M:%S) -- Alignment pipeline complete..." | tee -a $LOG

