#!/usr/bin/env Rscript
# Date:      8/20/2024  
# File:      utilities.R
# Author(s): Kasonde Chewe
# Description: Script includes custom functions and params for generation o


# Function to check if a package is installed, and if not, install it
install_pkg <- function(pkg, source = "CRAN") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (source == "Bioconductor") {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}


# Function to generate random 8-character alphanumeric strings
generate_uid_letter <- function(n) {
  replicate(n, paste0(sample(c(0:9, letters, LETTERS), 8, replace = TRUE), collapse = ""))
}


# Function to standardize and categorize Final Pathology
clean_pathology <- function(pathology) {
  pathology <- str_to_upper(pathology)  # Convert to uppercase
  
  # Reduce categories by pattern matching
  pathology <- case_when(
    str_detect(pathology, "PANCREATIC ADENOCARCINOMA|PDA") ~ "PANCREATIC ADENOCARCINOMA",
    str_detect(pathology, "DUCTAL ADENOCARCINOMA|DUCTAL") ~ "DUCTAL ADENOCARCINOMA",
    str_detect(pathology, "INVASIVE") ~ "INVASIVE ADENOCARCINOMA",
    str_detect(pathology, "METASTATIC") ~ "METASTATIC ADENOCARCINOMA",
    str_detect(pathology, "ADENOSQUAMOUS") ~ "ADENOSQUAMOUS CARCINOMA",
    str_detect(pathology, "POORLY DIFFERENTIATED") ~ "POORLY DIFFERENTIATED ADENOCARCINOMA",
    str_detect(pathology, "MODERATELY DIFFERENTIATED") ~ "MODERATELY DIFFERENTIATED ADENOCARCINOMA",
    str_detect(pathology, "WELL DIFFERENTIATED") ~ "WELL DIFFERENTIATED ADENOCARCINOMA",
    str_detect(pathology, "IPMN") ~ "INTRADUCTAL PAPILLARY MUCINOUS NEOPLASM (IPMN)",
    TRUE ~ pathology  # Keep original value if no pattern matched
  )
  
  return(pathology)
}


# Function to standardize and categorize Surgery
clean_surgery <- function(surgery) {
  surgery <- str_to_upper(surgery)  # Convert to uppercase
  
  # Reduce categories by pattern matching
  surgery <- case_when(
    str_detect(surgery, "WHIPPLE|PANCREATICODUODENECTOMY") ~ "PANCREATICODUODENECTOMY (WHIPPLE)",
    str_detect(surgery, "DISTAL PANCREATECTOMY|DISTAL/SPLEEN") ~ "DISTAL PANCREATECTOMY",
    str_detect(surgery, "TOTAL PANCREATECTOMY") ~ "TOTAL PANCREATECTOMY",
    str_detect(surgery, "EXPLORATORY LAPAROTOMY") ~ "EXPLORATORY LAPAROTOMY",
    str_detect(surgery, "FINE NEEDLE ASPIRATION|EUS-GUIDED BIOPSY") ~ "FINE NEEDLE ASPIRATION",
    str_detect(surgery, "LIVER BIOPSY") ~ "LIVER BIOPSY",
    str_detect(surgery, "ABORTED WHIPPLE") ~ "ABORTED WHIPPLE",
    str_detect(surgery, "OTHER|N/A|NA") ~ "OTHER/UNKNOWN",
    TRUE ~ surgery  # Keep original value if no pattern matched
  )
  
  return(surgery)
}


# Function to generate random 8-digit numeric strings
generate_uid_number <- function(n) {
  replicate(n, paste0(sample(0:9, 8, replace = TRUE), collapse = ""))
}