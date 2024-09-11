#!/usr/bin/env Rscript

# Author(s): Kasonde Chewe 
# Date: 09/10/2024
# Description: Script searches data repository and creates a manifest for WTS RNA files.

# Package Manager Function 
packages <- c("XML", "writexl")
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

install_and_load(packages)

# Set full_path to TRUE to include full file paths, or FALSE to include only filenames
full_path <- TRUE

# Set working directory
setwd("/oceanus/collab/InternalJeff/users/hxd052/data/JMP_PAC")

# Helper function to extract XML node values
get_node_value_by_name <- function(rootnode, parent_node, child_node) {
    parent <- rootnode[[parent_node]]
    if (!is.null(parent)) {
        child <- parent[[child_node]]
        if (!is.null(child)) {
            return(xmlValue(child))
        }
    }
    return(NA)
}

wts_df <- data.frame(JMPPAC_ID = character(),
                     CARIS_ID = character(),
                     WTS_RNA_BAM = character(),
                     WTS_RNA_FASTQ_R1 = character(),
                     WTS_RNA_FASTQ_R2 = character(),
                     WTS_RNA_TSV = character(),
                     WTS_REPORT_XML = character(),
                     stringsAsFactors = FALSE)

file_mode <- ifelse(full_path, TRUE, FALSE)

bam_files <- list.files(path = "WTS/BAM", pattern = "^RNA.*\\.bam$", full.names = file_mode)
fastq_r1_files <- list.files(path = "WTS/FASTQ", pattern = "^RNA.*.R1.*\\.fastq.gz$", full.names = file_mode)
fastq_r2_files <- list.files(path = "WTS/FASTQ", pattern = "^RNA.*.R2.*\\.fastq.gz$", full.names = file_mode)
tsv_files <- list.files(path = "WTS/exp", pattern = "^RNA.*\\.tsv$", full.names = file_mode)
xml_files <- list.files(path = "WES/reports/xml", pattern = "^TN.*\\.xml$", full.names = file_mode)

extract_caris_id_from_rna <- function(filename) {
    file_info <- strsplit(basename(filename), "-|_")[[1]]
    return(paste0(file_info[2], "-", file_info[3]))
}


for (bam_file in bam_files) {
    bam_filename <- ifelse(full_path, bam_file, basename(bam_file))
    CARIS_ID <- extract_caris_id_from_rna(bam_filename)
    wts_df <- rbind(wts_df, data.frame(JMPPAC_ID = "", CARIS_ID = CARIS_ID, WTS_RNA_BAM = bam_filename, 
                                       WTS_RNA_FASTQ_R1 = "", WTS_RNA_FASTQ_R2 = "", WTS_RNA_TSV = "", 
                                       WTS_REPORT_XML = "", stringsAsFactors = FALSE))
}

for (fastq_file in fastq_r1_files) {
    fastq_filename <- ifelse(full_path, fastq_file, basename(fastq_file))
    CARIS_ID <- extract_caris_id_from_rna(fastq_filename)
    wts_df$WTS_RNA_FASTQ_R1[wts_df$CARIS_ID == CARIS_ID] <- fastq_filename
}

for (fastq_file in fastq_r2_files) {
    fastq_filename <- ifelse(full_path, fastq_file, basename(fastq_file))
    CARIS_ID <- extract_caris_id_from_rna(fastq_filename)
    wts_df$WTS_RNA_FASTQ_R2[wts_df$CARIS_ID == CARIS_ID] <- fastq_filename
}

for (tsv_file in tsv_files) {
    tsv_filename <- ifelse(full_path, tsv_file, basename(tsv_file))
    CARIS_ID <- extract_caris_id_from_rna(tsv_filename)
    wts_df$WTS_RNA_TSV[wts_df$CARIS_ID == CARIS_ID] <- tsv_filename
}

for (xml_file in xml_files) {
    res <- xmlParse(file = xml_file)
    rootnode <- xmlRoot(res)
    caris_id_value <- get_node_value_by_name(rootnode, "testDetails", "labReportID")
    jmp_pac_num <- get_node_value_by_name(rootnode, "patientInformation", "firstName")
    jmp_pac_id <- paste0("JMP_Pac", jmp_pac_num)
    matching_rows <- grep(caris_id_value, wts_df$CARIS_ID)
    if (length(matching_rows) > 0) {
        wts_df$JMPPAC_ID[matching_rows] <- jmp_pac_id
        wts_df$WTS_REPORT_XML[matching_rows] <- ifelse(full_path, xml_file, basename(xml_file))
    }
}

# Print the final dataframe
print("Final dataframe:")
print(head(wts_df, 130))
print(dim(wts_df))

# Save the dataframe to an Excel or CSV file
write_xlsx(wts_df, "/oceanus/collab/InternalJeff/users/kxc732/Projects/../../hxd052/data/JMP_PAC/datasets/manifests/WTS_manifest.xlsx")
write.csv(wts_df, "/oceanus/collab/InternalJeff/users/kxc732/Projects/../../hxd052/data/JMP_PAC/datasets/manifests/WTS_manifest.csv")
