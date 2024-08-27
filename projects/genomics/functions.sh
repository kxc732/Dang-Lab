#!/bin/bash 

# Author:      Kasonde Chewe 
# File:        functions.sh
# Date:        8/16/2024
# Description: Custom functions for running neat NGS pipeline and report generation
# Note: analysis files should be isolated in per directory requirements 

# Updates 

echo "Loaded helper functions from helpers.sh"

# Function to extract sample names and file metadata, print a table, and/or save to CSV
extract_and_summarize_metadata() {
    local files_array=()
    local print_table="FALSE"
    local save_csv="FALSE"
    local output_csv=""

    # Parse command-line arguments
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            -f|--files_array) files_array=("${!2}"); shift ;;
            -p|--print_table) print_table="$2"; shift ;;
            -s|--save_csv) save_csv="TRUE"; output_csv="$2"; shift ;;
            -h|--help) echo "Usage: extract_and_summarize_metadata -f <files_array[@]> [-p <TRUE|FALSE>] [-s <output_csv_path>]"; return 0 ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if files_array is not empty
    if [ "${#files_array[@]}" -ne 0 ]; then
        log_message "Total fastq files found: ${#files_array[@]} !!!"
        log_message "Extracting sample names and file metadata"

        # Initialize the CSV file if save_csv is TRUE
        if [ "$save_csv" = "TRUE" ]; then
            if [ -z "$output_csv" ]; then
                echo "[ERROR] Output file path must be provided when save_csv is TRUE"
                return 1
            fi
            echo "Sample,FileSize,FilePath" > "$output_csv"
        fi

        # Initialize the summary table string
        local summary_table="Sample | FileSize | FilePath\n"

        # Loop through each file and gather the required information
        for fq_file in "${files_array[@]}"; do
            local sample=$(basename "$fq_file" | grep -oP '^[^\.]+')
            local filesize=$(du -sh "$fq_file" | cut -f1)
            local filepath="$fq_file"

            # Append the information to the summary_table string
            summary_table+="$sample | $filesize | $filepath\n"

            # Write the information to the CSV file if save_csv is TRUE
            if [ "$save_csv" = "TRUE" ]; then
                echo "$sample,$filesize,$filepath" >> "$output_csv"
            fi
        done

        # Print the summary table if print_table is TRUE
        if [ "$print_table" = "TRUE" ]; then
            printTable "|" "$summary_table"
        fi

        # Log the CSV generation if save_csv is TRUE
        if [ "$save_csv" = "TRUE" ]; then
            log_message "CSV file generated: $output_csv"
        fi
    else
        echo -e "\n -- No files found -- \n"
    fi
}



# Function to print formatted tables
printTable()
{
    local -r delimiter="${1}"
    local -r data="$(removeEmptyLines "${2}")"

    if [[ "${delimiter}" != '' && "$(isEmptyString "${data}")" = 'false' ]]
    then
        local -r numberOfLines="$(wc -l <<< "${data}")"

        if [[ "${numberOfLines}" -gt '0' ]]
        then
            local table=''
            local i=1

            for ((i = 1; i <= "${numberOfLines}"; i = i + 1))
            do
                local line=''
                line="$(sed "${i}q;d" <<< "${data}")"

                local numberOfColumns='0'
                numberOfColumns="$(awk -F "${delimiter}" '{print NF}' <<< "${line}")"

                # Add Line Delimiter
                if [[ "${i}" -eq '1' ]]
                then
                    table="${table}$(printf '%s#+' "$(repeatString '#+' "${numberOfColumns}")")"
                fi

                # Add Header Or Body
                table="${table}\n"

                local j=1
                for ((j = 1; j <= "${numberOfColumns}"; j = j + 1))
                do
                    table="${table}$(printf '#| %s' "$(cut -d "${delimiter}" -f "${j}" <<< "${line}")")"
                done

                table="${table}#|\n"

                # Add Line Delimiter
                if [[ "${i}" -eq '1' ]] || [[ "${numberOfLines}" -gt '1' && "${i}" -eq "${numberOfLines}" ]]
                then
                    table="${table}$(printf '%s#+' "$(repeatString '#+' "${numberOfColumns}")")"
                fi
            done

            if [[ "$(isEmptyString "${table}")" = 'false' ]]
            then
                echo -e "${table}" | column -s '#' -t | awk '/^\+/{gsub(" ", "-", $0)}1'
            fi
        fi
    fi
}

function removeEmptyLines()
{
    local -r content="${1}"
    echo -e "${content}" | sed '/^\s*$/d'
}

function repeatString()
{
    local -r string="${1}"
    local -r numberToRepeat="${2}"

    if [[ "${string}" != '' && "${numberToRepeat}" =~ ^[1-9][0-9]*$ ]]
    then
        local -r result="$(printf "%${numberToRepeat}s")"
        echo -e "${result// /${string}}"
    fi
}

function isEmptyString()
{
    local -r string="${1}"
    if [[ "$(trimString "${string}")" = '' ]]
    then
        echo 'true' && return 0
    fi
    echo 'false' && return 1
}

function trimString()
{
    local -r string="${1}"
    sed 's,^[[:blank:]]*,,' <<< "${string}" | sed 's,[[:blank:]]*$,,'
}

# Function to log messages with timestamp
log_message() {
    local message="$1"
    if [ $? -eq 0 ]; then
        echo "[LOG] -- $(date +%H:%M:%S) -- : $message"
    else
        echo "[ERROR] -- $(date +%H:%M:%S) -- : An error occurred on line $LINENO"
    fi
}


# Function to run FASTQC on an array of files
run_fastqc() {
    local -n files_array=$1
    local output_dir=$2

    mkdir -p $output_dir
    
    for fq_file in "${files_array[@]}"; do
        echo "[DEBUG] Running FASTQC on: $(basename $fq_file)"
        fastqc --threads 40 --svg -o "$output_dir" "$fq_file"
    done  
    log_message "FASTQC analysis completed. Results stored in $output_dir"
}



# Function to test if R1 and R2 files are correctly matched
test_paired_files() {
    files_array=("$@")
    matched=true

    for fq_file in "${files_array[@]}"; do
        # Check if the file is an R1 file
        if [[ "$fq_file" == *".r1"* ]]; then
            r1_file="$fq_file"
            r2_file="${fq_file/.r1/.r2}"
            
            # Check if the corresponding R2 file exists
            if [[ ! " ${files_array[*]} " =~ " $r2_file " ]]; then
                echo "[ERROR] Matching R2 file not found for $r1_file"
                matched=false
            fi
        fi
    done

    if [ "$matched" = true ]; then
        echo "[INFO] All files are correctly paired."
    else
        echo "[INFO] Some files are not correctly paired."
    fi
}


# Function to run Cutadapt on an array of files
run_cutadapt() {
    local files_array=("$@")
    local output_dir="${files_array[-1]}"
    unset files_array[-1]  # Remove the last element, which is the output directory
    
    local adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    local adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Loop through each file in the array and run Cutadapt
    for fq_file in "${files_array[@]}"; do
        # Check if the file is an R1 file
        if [[ "$fq_file" == *".r1"* ]]; then
            r1_file="$fq_file"
            r2_file="${fq_file/.r1/.r2}"
            
            # Ensure the R2 file exists
            if [[ ! " ${files_array[*]} " =~ " $r2_file " ]]; then
                echo "[ERROR] Matching R2 file not found for $r1_file"
                continue
            fi
            
            # Get the base name of the file for naming output
            base_name=$(basename "$r1_file" .r1.fastq.gz)
            
            # Run Cutadapt
            cutadapt \
                -a "$adapter1" \
                -A "$adapter2" \
                -q 20 \
                -m 30 \
                -j 0 \
                -o "$output_dir/${base_name}_trimmed_R1.fastq.gz" \
                -p "$output_dir/${base_name}_trimmed_R2.fastq.gz" \
                "$r1_file" \
                "$r2_file"
        fi
    done

    log_message "Cutadapt trimming completed. Trimmed files stored in $output_dir"
}


# Function to find the best matching R2 file for a given R1 file
find_matching_r2() {
    local r1_file=$1
    local -n all_files_array=$2  # Reference to the array of all files

    # Extract the base name from the R1 file, replacing _R1 with _R2
    local r2_file_base=$(basename "$r1_file" | sed 's/_R1/_R2/')

    # Loop through all files to find the matching R2 file
    for fq_file in "${all_files_array[@]}"; do
        if [[ "$(basename "$fq_file")" == "$r2_file_base" ]]; then
            echo "$fq_file"
            return 0
        fi
    done

    echo "[ERROR] No match found for $r1_file"
    return 1
}



# Function to run fastp on an array of files
run_fastp() {
    local -n files_array=$1  # Reference to the array of files
    local output_dir=$2       # Output directory

    mkdir -p $output_dir

    # First, find all R1 files and their matching R2 files
    for fq_file in "${files_array[@]}"; do
        log_message "Processing file: $fq_file"
        
        # Convert to lowercase and check if the file has the .R1.fastq.gz suffix
        lower_fq_file=$(echo "$fq_file" | tr '[:upper:]' '[:lower:]')
        if [[ "$lower_fq_file" == *".r1.fastq.gz" ]]; then
            r1_file="$fq_file"
            log_message "Identified R1 file: $r1_file"
            
            # Capture the output of the find_matching_r2 function in the r2_file variable
            r2_file=$(find_matching_r2 "$r1_file" files_array)
            
            # Ensure the R2 file exists
            if [[ -z "$r2_file" ]]; then
                log_message "[ERROR] Matching R2 file not found for $r1_file"
                continue
            fi

            log_message "Identified matching R2 file: $r2_file"

            # Get the base name of the file for naming output
            base_name=$(basename "$r1_file" | sed 's/\.r1\.fastq\.gz//I')
            log_message "Base name for output files: $base_name"

            # Run fastp with 16 threads
            log_message "Running fastp on $r1_file and $r2_file"
            fastp \
                -i "$r1_file" \
                -I "$r2_file" \
                -o "$output_dir/${base_name}_trimmed_R1.fastq.gz" \
                -O "$output_dir/${base_name}_trimmed_R2.fastq.gz" \
                --html "$output_dir/${base_name}_fastp_report.html" \
                --json "$output_dir/${base_name}_fastp_report.json" \
                --report_title "$base_name fastp report" \
                --qualified_quality_phred 20 \
                --length_required 30 \
                --thread 16
        fi
    done

    log_message "fastp trimming completed. Trimmed files and reports stored in $output_dir"
}



# Function to create summary reports from fastp JSON output
summarize_fastp_reports() {
    local json_files=()
    local output_csv=""
    local print_table="FALSE"

    # Parse the input arguments
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --output_csv) output_csv="$2"; shift ;;
            --print_table) print_table="$2"; shift ;;
            -f|--files) shift ;;  # Skip files flag
            *) json_files+=("$1") ;;  # Add JSON files to the array
        esac
        shift
    done

    if [ -z "$output_csv" ]; then
        log_message "[ERROR] Output CSV path must be specified with --output_csv"
        return 1
    fi

    # Initialize the CSV with headers
    echo "sample_name,fastp_version,sequencing,before_filtering_total_reads,before_filtering_total_bases,before_filtering_q20_bases,before_filtering_q30_bases,before_filtering_q20_rate,before_filtering_q30_rate,before_filtering_read1_mean_length,before_filtering_read2_mean_length,before_filtering_gc_content,after_filtering_total_reads,after_filtering_total_bases,after_filtering_q20_bases,after_filtering_q30_bases,after_filtering_q20_rate,after_filtering_q30_rate,after_filtering_read1_mean_length,after_filtering_read2_mean_length,after_filtering_gc_content" > "$output_csv"
    log_message "Initialized CSV file with headers: $output_csv"

    # Initialize the summary table string
    local summary_table="sample_name | fastp_version | sequencing | before_filtering_total_reads | before_filtering_total_bases | before_filtering_q20_bases | before_filtering_q30_bases | before_filtering_q20_rate | before_filtering_q30_rate | before_filtering_read1_mean_length | before_filtering_read2_mean_length | before_filtering_gc_content | after_filtering_total_reads | after_filtering_total_bases | after_filtering_q20_bases | after_filtering_q30_bases | after_filtering_q20_rate | after_filtering_q30_rate | after_filtering_read1_mean_length | after_filtering_read2_mean_length | after_filtering_gc_content\n"

    # Loop through each JSON file
    for json_file in "${json_files[@]}"; do
        # Extract the sample name from the filename
        sample_name=$(basename "$json_file" | sed 's/_fastp_report.json//')

        # Extract the data using jq
        row=$(jq -r --arg sample_name "$sample_name" '
        [
            $sample_name,
            .summary.fastp_version, 
            .summary.sequencing, 
            .summary.before_filtering.total_reads, 
            .summary.before_filtering.total_bases, 
            .summary.before_filtering.q20_bases, 
            .summary.before_filtering.q30_bases, 
            .summary.before_filtering.q20_rate, 
            .summary.before_filtering.q30_rate, 
            .summary.before_filtering.read1_mean_length, 
            .summary.before_filtering.read2_mean_length, 
            .summary.before_filtering.gc_content,
            .summary.after_filtering.total_reads, 
            .summary.after_filtering.total_bases, 
            .summary.after_filtering.q20_bases, 
            .summary.after_filtering.q30_bases, 
            .summary.after_filtering.q20_rate, 
            .summary.after_filtering.q30_rate, 
            .summary.after_filtering.read1_mean_length, 
            .summary.after_filtering.read2_mean_length, 
            .summary.after_filtering.gc_content
        ]
        | @csv' "$json_file")

        # Append the row to the CSV
        echo "$row" >> "$output_csv"

        # Append the row to the summary table string
        summary_table+=$(echo "$row" | tr ',' ' | ')$'\n'
        
        log_message "Processed ---- "${json_file##*/}" "
    done

    log_message "Summary created in $output_csv"

    # Print the summary table if requested
    if [ "$print_table" = "TRUE" ]; then
        log_message "Printing the summary table to the terminal"
        printTable "|" "$summary_table"
    fi
}




# Function to perform peak calling using MACS2
run_peakcalling() {
    # to-do: send each instance of a command to slurm to.. that is add param if slurm 
    # to-do: add a default threads param

    local bam_files=()  # Array of BAM files
    local output_dir=""  # Output directory
    local genome_size="hs"  # Default genome size (human)
    local pvalue=""  # Default p-value threshold for peak calling
    local qvalue=""  # Default q-value (FDR)
    local input_normalization="FALSE"  # Default to not use input normalization
    local binsize=1  # Default bin size

    # Parse command-line arguments
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            -b|--bam_files) bam_files=("${!2}"); shift ;;
            -o|--output_dir) output_dir="$2"; shift ;;
            -g|--genome_size) genome_size="$2"; shift ;;
            -p|--pvalue) pvalue="$2"; shift ;;
            -q|--qvalue) qvalue="$2"; shift ;;
            --input-normalization) input_normalization="$2"; shift ;;
            --binsize) binsize="$2"; shift ;;
            -h|--help) 
                echo "Usage: run_peakcalling -b <bam_files[@]> -o <output_dir> [-g <genome_size>] [-p <pvalue>] [-q <qvalue>] [--input-normalization <TRUE|FALSE>] [--binsize <size>]"
                return 0 
                ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required parameters are provided
    if [ -z "$output_dir" ] || [ "${#bam_files[@]}" -eq 0 ]; then
        echo "[ERROR] Both --bam_files and --output_dir must be provided."
        return 1
    fi

    # Ensure only one of qvalue or pvalue is set
    if [ -n "$qvalue" ] && [ -n "$pvalue" ]; then
        echo "[ERROR] Only one of --qvalue or --pvalue should be provided."
        return 1
    fi

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Loop through the array of BAM files and run MACS2 on each
    for bam_file in "${bam_files[@]}"; do
        # Get the base name of the BAM file for naming outputs
        local sample_name=$(basename "$bam_file" .final.bam)

        # Skip input control files when input normalization is enabled
        if [ "$input_normalization" = "TRUE" ] && [[ "$sample_name" == *input_trimmed ]]; then
            # echo "[DEBUG] Skipping input file $sample_name since input normalization is enabled."
            continue
        fi

        # Initialize control BAM file
        local control_bam=""

        # If input normalization is enabled, find the corresponding control BAM
        if [ "$input_normalization" = "TRUE" ]; then
            local input_name="${sample_name/_trimmed/-input_trimmed}.final.bam"
            # echo "[DEBUG] Looking for input BAM: $input_name"

            # Attempt to find the control BAM in the provided bam_files
            for file in "${bam_files[@]}"; do
                if [[ "$file" == *"$input_name" ]]; then
                    control_bam="$file"
                    break
                fi
            done

            if [ -z "$control_bam" ]; then
                echo "[ERROR] No matching input control BAM file found for $sample_name"
                continue
            fi
        fi

        # Construct the MACS2 command
        local macs2_cmd="macs2 callpeak -t \"$bam_file\" -f BAM -g \"$genome_size\"  -f \"BAMPE\" -n \"$sample_name\" --outdir \"$output_dir\" --nomodel --extsize \"$binsize\" --keep-dup all"

        # Add either pvalue or qvalue to the MACS2 command
        if [ -n "$pvalue" ]; then
            macs2_cmd+=" --pvalue \"$pvalue\""
        elif [ -n "$qvalue" ]; then
            macs2_cmd+=" --qvalue \"$qvalue\""
        fi

        # Add the control BAM to the command if it exists
        if [ -n "$control_bam" ]; then
            macs2_cmd+=" -c \"$control_bam\""
        fi

        # Run the MACS2 command
        eval "$macs2_cmd"
        # Log completion message
        log_message "MACS2 peak calling completed for $bam_file. Output stored in $output_dir/$sample_name"
    done
}


# Function to convert bams to bw
run_bam_to_bw () {
    local bams_array=()
    local output_dir=""
    local binsize=1       # Default bin size 
    local threads=40      # Default use 64 threads
    local effective_genome_size=2913022398  
    local normalize_using="RPGC"  # Default normalization method

    # Parse command-line arguments
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            -b|--bams_array) bams_array=("${!2}"); shift ;;
            -o|--output_dir) output_dir="$2"; shift ;;
            -t|--threads) threads="$2"; shift 2 ;;
            --binsize) binsize="$2"; shift 2 ;;
            --effectiveGenomeSize) effective_genome_size="$2"; shift 2 ;;
            --normalizeUsing) normalize_using="$2"; shift 2 ;;
            -h|--help) 
                echo "Usage: run_bam_to_bw -f <bams_array[@]> -o <output_dir> -r <reference_genome> -b <blacklist> [-t <TRUE|FALSE>]"; 
                return 0 ;;
            *) echo "Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if required arguments are provided
    if [[ ${#bams_array[@]} -eq 0 || -z "$output_dir" ]]; then
        echo "[ERROR] Both --bam_files and --output_dir must be provided."
        return 1
    fi

    # Loop through each BAM file and run bamCoverage
    for bam_file in "${bams_array[@]}"; do
        base_name=$(basename "$bam_file" .bam)
        output_bw="${output_dir}/${base_name}.bw"

        # Execute bamCoverage command on slurm 
        bamCoverage -b "$bam_file" -o "$output_bw" --binSize "$binsize" -p "$threads" --effectiveGenomeSize "$effective_genome_size" --normalizeUsing "$normalize_using"
    done
}

# Function to run BWA alignment on BAM files
run_bwa_alignment() {
    local files_array=()
    local output_dir=""
    local reference_genome=""
    local blacklist=""
    local test_mode="FALSE"

    # Parse command-line arguments
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            -f|--files_array) files_array=("${!2}"); shift ;;
            -o|--output_dir) output_dir="$2"; shift ;;
            -r|--reference_genome) reference_genome="$2"; shift ;;
            -b|--blacklist) blacklist="$2"; shift ;;
            -t|--test_mode) test_mode="$2"; shift ;;
            -h|--help) echo "Usage: run_bwa_alignment -f <files_array[@]> -o <output_dir> -r <reference_genome> -b <blacklist> [-t <TRUE|FALSE>]"; return 0 ;;
            *) echo "[ERROR] Unknown parameter passed: $1"; return 1 ;;
        esac
        shift
    done

    # Check if all required parameters are provided
    if [ -z "$output_dir" ] || [ -z "$reference_genome" ] || [ -z "$blacklist" ] || [ "${#files_array[@]}" -eq 0 ]; then
        echo "[ERROR] All parameters (--files_array, --output_dir, --reference_genome, --blacklist) must be provided."
        return 1
    fi

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Counter to check if we're in test mode and only processing the first pair
    local processed_pairs=0

    # Loop through the array of files
    for fq_file in "${files_array[@]}"; do
        
        # Check if the file has the _R1.fastq.gz suffixA
        if [[ "$fq_file" == *_R1.fastq.gz ]]; then
            r1_file="$fq_file"            
            # Use find_matching_r2 to find the corresponding R2 file
            r2_file=$(find_matching_r2 "$r1_file" files_array)
            
            # Ensure the corresponding R2 file exists
            if [[ -z "$r2_file" ]]; then
                echo "[ERROR] Matching R2 file not found for $r1_file"
                continue
            fi

            # Get the base name of the file for naming output
            # To-do: automate to remove suffix... perhaps dictionary of alternative suffixes
            base_name=$(basename "$r1_file" | sed 's/_trimmed_R1.fastq.gz//')
            log_message "Initializing alignment for: $base_name"

            # Create sample output directory
            local sample_output_dir="$output_dir/$base_name"
            mkdir -p "$sample_output_dir"

            # Run BWA MEM
            log_message "Running BWA MEM on $r1_file and $r2_file"
            bwa mem -t 60 -c 1 "$reference_genome" "$r1_file" "$r2_file" > "$sample_output_dir/$base_name.sam"
            log_message "BWA alignment completed for $r1_file and $r2_file. Output: $sample_output_dir/$base_name.sam"

            # Convert SAM to sorted BAM
            log_message "Converting SAM to sorted BAM for $base_name"
            samtools view -b -@ 60 "$sample_output_dir/$base_name.sam" | samtools sort -@ 60 -o "$sample_output_dir/$base_name.sorted.bam" -O bam
            samtools index "$sample_output_dir/$base_name.sorted.bam"
            log_message "SAM to sorted BAM conversion completed for $base_name. Output: $sample_output_dir/$base_name.sorted.bam"

            # Filter out mitochondrial reads
            log_message "Filtering out mitochondrial reads for $base_name"
            samtools view -h "$sample_output_dir/$base_name.sorted.bam" | grep -v 'chrM' | samtools view -b -@ 60 > "$sample_output_dir/$base_name.sorted.chrM.bam"
            log_message "Filtered out mitochondrial reads for $base_name. Output: $sample_output_dir/$base_name.sorted.chrM.bam"

            # Mark duplicates
            log_message "Marking duplicates for $base_name"
            picard MarkDuplicates -I "$sample_output_dir/$base_name.sorted.chrM.bam" -M "$sample_output_dir/$base_name.metrics.txt" -O "$sample_output_dir/$base_name.sorted.chrM.dups.bam" --REMOVE_DUPLICATES true
            log_message "Marked duplicates for $base_name. Output: $sample_output_dir/$base_name.sorted.chrM.dups.bam"

            # Filter BAM file
            log_message "Filtering BAM file for $base_name"
            samtools view -b -@ 60 -F 3076 -q 30 -o "$sample_output_dir/$base_name.filtered.bam" "$sample_output_dir/$base_name.sorted.chrM.dups.bam"
            samtools index "$sample_output_dir/$base_name.filtered.bam"
            log_message "Filtered BAM file for $base_name. Output: $sample_output_dir/$base_name.filtered.bam"

            # Remove blacklisted regions
            log_message "Removing blacklisted regions for $base_name"
            bedtools intersect -v -abam "$sample_output_dir/$base_name.filtered.bam" -b "$blacklist" > "$sample_output_dir/$base_name.final.bam"
            samtools index "$sample_output_dir/$base_name.final.bam"
            log_message "Removed blacklisted regions for $base_name. Output: $sample_output_dir/$base_name.final.bam"

            # Compress the SAM file
            log_message "Compressing SAM file for $base_name"
            pigz "$sample_output_dir/$base_name.sam"
            log_message "Compressed SAM file for $base_name. Output: $sample_output_dir/$base_name.sam.gz"

            # Increment the processed pairs counter
            processed_pairs=$((processed_pairs + 1))
            log_message "Processed pairs count: $processed_pairs"

            # If test mode is enabled and we've processed the first pair, exit the loop
            if [ "$test_mode" == "TRUE" ] && [ "$processed_pairs" -ge 1 ]; then
                log_message "Test mode is enabled. Only the first pair of files was processed."
                break
            fi
        fi
    done
    log_message "BWA alignment process completed."
}





# Function to normalize BigWig files by matching IP and INPUT files using bigwigCompare
run_bw_normalization() {
    local ip_files=("${!1}")  # Array of IP files
    local input_files=("${!2}")  # Array of INPUT files
    local output_dir="$3"  # Output directory

    # Ensure the output directory exists
    mkdir -p "$output_dir"

    # Loop through each IP file
    for ip_file in "${ip_files[@]}"; do
        local ip_basename=$(basename "$ip_file")
        local ip_prefix=${ip_basename%-*}  # Strip everything after the last '-'

        # Initialize variables to track the best match
        local best_match=""
        local best_prefix_match_length=0

        # Loop through each INPUT file to find the best match
        for input_file in "${input_files[@]}"; do
            local input_basename=$(basename "$input_file")
            local input_prefix=${input_basename%-*}  # Strip everything after the last '-'

            # Check how much of the prefixes match
            local common_prefix_length=${#ip_prefix}
            while [[ "$common_prefix_length" -gt 0 && "${ip_prefix:0:$common_prefix_length}" != "${input_prefix:0:$common_prefix_length}" ]]; do
                ((common_prefix_length--))
            done

            # If this is the best match so far, update best_match
            if [[ $common_prefix_length -gt $best_prefix_match_length ]]; then
                best_match="$input_file"
                best_prefix_match_length=$common_prefix_length
            fi
        done

        # Resolve final output filename as <prefix>_norm.bw
        local final_output="$output_dir/${ip_prefix}_norm.bw"

        # Perform normalization using bigwigCompare
        echo "Normalizing $ip_file with $best_match"
        bigwigCompare -b1 "$ip_file" -b2 "$best_match" -o "$final_output" -of bigwig -bs 1 -p 77 --operation subtract

        echo "Output saved as $final_output"
    done
}


# Function to get sample names by stripping the suffix with debug prints
get_samplenames() {
    local files=("$@")  # Array of files passed as the first argument
    local suffix=$1     # Suffix to strip passed as the second argument
    shift               # Remove the first argument (suffix) from the array
    local stripped_names=()  # Array to store the stripped names

    for file in "${files[@]}"; do

        basename=$(basename "$file")  # Get the basename of the file
        stripped_name=${basename%"$suffix"}  # Strip the suffix from the basename
        stripped_names+=("$stripped_name")  # Add the stripped name to the array
    done
    echo "${stripped_names[@]}"  # Return the array of stripped names
}

# Unused function keep for archive purposes
# Function to match IP and INPUT files
match_ip_input_files() {
    local ip_files=("${!1}")  # Array of IP files
    local input_files=("${!2}")  # Array of INPUT files
    declare -A matched_pairs  # Associative array to store matched pairs

    # Loop through each IP file
    for ip_file in "${ip_files[@]}"; do
        local ip_basename=$(basename "$ip_file")
        local ip_prefix=${ip_basename%-*}  # Strip everything after the last '-'

        # Initialize variables to track the best match
        local best_match=""
        local best_prefix_match_length=0

        # Loop through each INPUT file
        for input_file in "${input_files[@]}"; do
            local input_basename=$(basename "$input_file")
            local input_prefix=${input_basename%-*}  # Strip everything after the last '-'

            # Check how much of the prefixes match
            local common_prefix_length=${#ip_prefix}
            while [[ "$common_prefix_length" -gt 0 && "${ip_prefix:0:$common_prefix_length}" != "${input_prefix:0:$common_prefix_length}" ]]; do
                ((common_prefix_length--))
            done

            # If this is the best match so far, update best_match
            if [[ $common_prefix_length -gt $best_prefix_match_length ]]; then
                best_match="$input_file"
                best_prefix_match_length=$common_prefix_length
            fi
        done

        # Store the best match
        matched_pairs["$ip_file"]="$best_match"
    done

    # Output the matched pairs
    for ip_file in "${!matched_pairs[@]}"; do
        echo "IP file: $ip_file"
        echo "Matched INPUT file: ${matched_pairs[$ip_file]}"
        echo
    done
}

generate_bam_stats() {
    local bam_files=("$@")  # Array of BAM files
    local output_dir="${bam_files[-1]}"  # The last element is the output directory
    unset bam_files[-1]  # Remove the last element from the array

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Initialize the CSV file with headers
    local output_csv="$output_dir/bam_alignment_stats.csv"
    echo "sample_name,total_reads,mapped_reads,mapped_rate,average_length,insert_size_mean,insert_size_std_dev" > "$output_csv"
    log_message "Initialized CSV file with headers: $output_csv"

    # # Loop through each BAM file and generate stats
    # for bam_file in "${bam_files[@]}"; do
    #     # Extract the sample name from the BAM file name
    #     local sample_name=$(basename "$bam_file" .bam)

    #     # Run samtools stats and capture the output
    #     local stats=$(samtools stats "$bam_file")

    #     # Extract key metrics using grep and awk
    #     local total_reads=$(echo "$stats" | grep "^SN" | grep "raw total sequences:" | awk '{print $4}')
    #     local mapped_reads=$(echo "$stats" | grep "^SN" | grep "reads mapped:" | awk '{print $4}')
    #     local mapped_rate=$(echo "$stats" | grep "^SN" | grep "reads mapped:" | awk '{print $5}' | tr -d '()')
    #     local average_length=$(echo "$stats" | grep "^SN" | grep "average length:" | awk '{print $4}')
    #     local insert_size_mean=$(echo "$stats" | grep "^IS" | grep "insert size average:" | awk '{print $4}')
    #     local insert_size_std_dev=$(echo "$stats" | grep "^IS" | grep "insert size standard deviation:" | awk '{print $5}')

    #     # Append the metrics to the CSV file
    #     echo "$sample_name,$total_reads,$mapped_reads,$mapped_rate,$average_length,$insert_size_mean,$insert_size_std_dev" >> "$output_csv"
    #     log_message "Processed BAM file: $bam_file"
    # done

    log_message "BAM stats collection completed. Combined stats CSV generated: $output_csv"
}
