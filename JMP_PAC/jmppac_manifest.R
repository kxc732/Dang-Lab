#!/usr/bin/env Rscript

# Author(s):
# Date: 09/10/2024
# Description: Script searches data repository an


# Package Manager Function 
packages <- c("XML", "VennDiagram", "writexl")
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

install_and_load(packages)


# Working directory to project files
setwd("/oceanus/collab/InternalJeff/users/hxd052/data/JMP_PAC")

# Extract Whole Exome Sequencing (WES) Manifest ---
# --------------------------------------------------
setwd("/oceanus/collab/InternalJeff/users/hxd052/data/JMP_PAC")

if (!require(XML)) {
    install.packages("XML")
    library(XML)
}

if (!require(writexl)) {
    install.packages("writexl")
    library(writexl)
}

files <- list.files(recursive = TRUE)
print(files)

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

df <- data.frame(JMPPAC_ID = character(),
                 CARIS_ID = character(),
                 WES_DNA_BAM = character(),
                 WES_DNA_FASTQ_R1 = character(),
                 WES_DNA_FASTQ_R2 = character(),
                 WES_DNA_VCF = character(),
                 WES_REPORT_PDF = character(),
                 WES_REPORT_XML = character(),
                 WES_REPORT_JSON = character(),
                 stringsAsFactors = FALSE)

print(df)

bam_files <- list.files(path = "WES/BAM", pattern = "^DNA.*\\.bam$", full.names = TRUE)
fastq_r1_files <- list.files(path = "WES/FASTQ", pattern = "^DNA.*.R1.*\\.fastq.gz$", full.names = TRUE)
fastq_r2_files <- list.files(path = "WES/FASTQ", pattern = "^DNA.*.R2.*\\.fastq.gz$", full.names = TRUE)
vcf_files <- list.files(path = "WES/vcf", pattern = "^DNA.*\\.vcf$", full.names = TRUE)
pdf_files <- list.files(path = "WES/reports/pdf", pattern = "^TN.*\\.pdf$", full.names = TRUE)
xml_files <- list.files(path = "WES/reports/xml", pattern = "^TN.*\\.xml$", full.names = TRUE)
json_files <- list.files(path = "WES/reports/json", pattern = "^TN.*\\.json$", full.names = TRUE)

extract_caris_id_from_tn <- function(filename) {
    file_info <- strsplit(filename, "-|_")[[1]]
    return(file_info[2])
}

for (bam_file in bam_files) {
    bam_filename <- basename(bam_file)
    file_info <- strsplit(bam_filename, "_|\\.")[[1]]
    CARIS_ID <- file_info[2]
    df <- rbind(df, data.frame(JMPPAC_ID = "", CARIS_ID = CARIS_ID, WES_DNA_BAM = bam_filename, WES_DNA_FASTQ_R1 = "", WES_DNA_FASTQ_R2 = "", WES_DNA_VCF = "", WES_REPORT_PDF = "", WES_REPORT_XML = "", stringsAsFactors = FALSE))
}

for (fastq_file in fastq_r1_files) {
    fastq_filename <- basename(fastq_file)
    file_info <- strsplit(fastq_filename, "_|\\.")[[1]]
    CARIS_ID <- file_info[2]
    df$WES_DNA_FASTQ_R1[df$CARIS_ID == CARIS_ID] <- fastq_filename
}

for (fastq_file in fastq_r2_files) {
    fastq_filename <- basename(fastq_file)
    file_info <- strsplit(fastq_filename, "_|\\.")[[1]]
    CARIS_ID <- file_info[2]
    df$WES_DNA_FASTQ_R2[df$CARIS_ID == CARIS_ID] <- fastq_filename
}

for (vcf_file in vcf_files) {
    vcf_filename <- basename(vcf_file)
    file_info <- strsplit(vcf_filename, "_|\\.")[[1]]
    CARIS_ID <- file_info[2]
    df$WES_DNA_VCF[df$CARIS_ID == CARIS_ID] <- vcf_filename
}

for (pdf_file in pdf_files) {
    pdf_filename <- basename(pdf_file)
    CARIS_ID <- extract_caris_id_from_tn(pdf_filename)
    matching_rows <- grep(CARIS_ID, df$CARIS_ID)
    if (length(matching_rows) > 0) {
        df$WES_REPORT_PDF[matching_rows] <- pdf_filename
    }
}

for (xml_file in xml_files) {
    xml_filename <- basename(xml_file)
    CARIS_ID <- extract_caris_id_from_tn(xml_filename)
    matching_rows <- grep(CARIS_ID, df$CARIS_ID)
    if (length(matching_rows) > 0) {
        df$WES_REPORT_XML[matching_rows] <- xml_filename
    }
}

for (json_file in json_files) {
    json_filename <- basename(json_file)
    CARIS_ID <- extract_caris_id_from_tn(json_filename)
    matching_rows <- grep(CARIS_ID, df$CARIS_ID)
    if (length(matching_rows) > 0) {
        df$WES_REPORT_JSON[matching_rows] <- json_filename
    }
}

head(df)

for (xml_file in xml_files) {
    res <- xmlParse(file = xml_file)
    rootnode <- xmlRoot(res)
    caris_id_value <- get_node_value_by_name(rootnode, "testDetails", "labReportID")
    jmp_pac_num <- get_node_value_by_name(rootnode, "patientInformation", "firstName")
    jmp_pac_id <- paste0("JMP_Pac", jmp_pac_num)
    matching_rows <- grep(caris_id_value, df$CARIS_ID)
    if (length(matching_rows) > 0) {
        df$JMPPAC_ID[matching_rows] <- jmp_pac_id
    }
}

#write_xlsx(df, "/oceanus/collab/InternalJeff/users/kxc732/Projects/../../hxd052/data/JMP_PAC/datasets/manifests/WES_manifest.xlsx")


# Print the final dataframe
print("Final dataframe:")
print(head(df))

