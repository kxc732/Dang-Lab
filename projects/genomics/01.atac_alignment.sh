#!/bin/bash

# FILE INFORMATION
#---------------------------
# Author:      Kasonde Chewe
# Pipeline:    ATAC-SEQ
# File:        01.atac_alignment.sh
# Date:        08/26/2024
# Description: Script for processing ATAC-SEQ data including FASTQ file processing, QC, trimming, and report generation.


# SLURM DIRECTIVES
#---------------------------
#SBATCH --job-name=ATACSEQ
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --nodes=5
#SBATCH --ntasks=40
#SBATCH --mem=128G

# ACTIVATE CONDA ENVIRONMENT
#----------------------------
eval "$(conda shell.bash hook)"
conda activate gro-seq

# SOURCE FUNCTIONS
# ---------------------------
source ~/Projects/Condensate_Paper/05_SMARC_FINAL/00_SCRIPTS_SMARC/utilities/bash/functions.sh

# VARIABLES AND PATHS
# ---------------------------
fastq_dir=~/Data_Repository/00ATAC-SEQ/raw # FASTQ directory
output_dir=~/Projects/Condensate_Paper/00_ATACSEQ_FINAL # this is working dir

# ENSURE DIRECTORIES EXIST 
# --------------------------
mkdir -p $output_dir/datasets

# Collect FASTQ files
fq_files=($(find "$fastq_dir" -type f -name "*.fastq.gz"))

# REPORT 1: FASTQ FILES PRE-PROCESSED SIZE AND PATHS
# --------------------------------------------------
output_csv1="$output_dir/datasets/00_atac_fastq_metadata.csv"
log_message "Generating metadata report for pre-processed FASTQ files."
extract_and_summarize_metadata -f fq_files[@] -p FALSE -s "$output_csv1"

# PRE-TRIM FASTQC 
# -----------------------------------------------------
fastqc_dir="$output_dir/fastqc"
mkdir -p $fastqc_dir mkdir -p $fastqc_dir/pre_trimmed
log_message "Running Pre-Trim FASTQC"
# run_fastqc fq_files "$fastqc_dir/pre_trimmed"

# TRIMMING
# -----------------------------------------------------
trimmed_dir="$output_dir/fastp"
mkdir -p $trimmed_dir
log_message "Running FASTP -- trimming FASTQ files"
run_fastp fq_files "$trimmed_dir"

# REPORT 2: TRIMMED FASTQ FILES SIZE AND PATHS
# -----------------------------------------------------
fq_files_trimmed=($(find "$trimmed_dir" -type f -name "*.fastq.gz*"))
output_csv2="$output_dir/datasets/01_atac_trimmed_fastq_metadata.csv"
log_message "Generating metadata report for trimmed FASTQ files."
extract_and_summarize_metadata -f fq_files_trimmed[@] -p TRUE -s "$output_csv2" 

# REPORT 3: TRIMMING STATISTICS FROM FASTP JSON
# -----------------------------------------------------
json_fq_trimmed=($(find "$trimmed_dir" -type f -name "*.json"))
output_csv3="$output_dir/datasets/02_atac_trimmed_fastq_statistics.csv"
log_message "Generating trimming statistics report from FASTP JSON files."
summarize_fastp_reports --output_csv "$output_csv3" --print_table "FALSE" "${json_fq_trimmed[@]}"

# POST-TRIM FASTQC
# -----------------------------------------------------
log_message "Running Post-Trim FASTQC"
run_fastqc fq_files_trimmed "$trimmed_dir"

# ALIGN READS TO REFERENCE GENOME HG38
# -------------------------------------------------------
reference_genome=~/genomes/hg38/index/bwa/hg38.fa
blacklist=~/projects-2024/scripts/utilities/blacklist_regions/RegionsRemove.bed
mapped_dir=$output_dir/mapped
mkdir -p $mapped_dir $mapped_dir/00bams


log_message "Running Alignment"
run_bwa_alignment -f fq_files_trimmed[@] -o $mapped_dir/00bams -r $reference_genome -b $blacklist 
log_message " ================= ATAC-SEQ alignment pipeline complete!!!!!!!!!! ================ "