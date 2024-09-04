#!/usr/bin/env Rscript
# Date: 09/04/2024
# Author(s): Kasonde Chewe
# Description: Script to deidentify JMP PAC dataset


# Libraries ----
library(ComplexHeatmap)
library(readxl)
library(lubridate)
library(tidyverse)
library(forcats)
library(dplyr)
library(grid)


# Variables and functions  ----
setwd("~/Projects/JMP_PAC/00_SCRIPTS_R/bin/R/final")
source("utilities.R")



# Clean JMP PAC MST LIST ----
# 1. Load dataset 
mst_list <- read_xlsx("~/Projects/JMP_PAC/datasets/pre/JMP_PaC_NAFLD_n502_Extended_Dataset.xlsx")


# 2. Drop Identifier Columns 
mst_list <- mst_list %>% dplyr::select(-NAME,-MRN, -`SPECIMEN ID`, -`SURGEON`,
                                       -`Recurrence Biopsy Specimen ID`, 
                                       -`Comments/Notes`)

# 3. Case relevant columns to date format in MM/DD/YYYY
mst_list <- mst_list %>%
  # a. Cast relevant columns to date format in MM/DD/YYYY
  mutate(
    DOB = format(as.Date(DOB), "%m/%d/%Y"),
    DOS = ifelse(
      grepl("^[0-9]+$", DOS) & !is.na(DOS),  # Check if DOS is numeric and not NA
      format(as.Date(as.numeric(DOS), origin = "1899-12-30"), "%m/%d/%Y"),
      NA  # Assign NA if DOS is not numeric
    ),
    `Date of Diagnosis` = ifelse(
      grepl("^[0-9]+$", `Date of Diagnosis`) & !is.na(`Date of Diagnosis`),  # Check if numeric and not NA
      format(as.Date(as.numeric(`Date of Diagnosis`), origin = "1899-12-30"), "%m/%d/%Y"),
      NA  # Assign NA if non-numeric
    )
  ) %>%
  # b. Correctly Populate Column Age at Diagnosis as integer
  mutate(
    `Age at Diagnosis` = ifelse(
      !is.na(DOS) & as.Date(DOS, "%m/%d/%Y") >= as.Date("01/01/2022", "%m/%d/%Y"),  # If DOS is valid and from 2022 and beyond
      as.integer(round(as.numeric(difftime(as.Date(DOS, "%m/%d/%Y"), as.Date(DOB, "%m/%d/%Y"), units = "days")) / 365.25)),
      ifelse(
        !is.na(`Date of Diagnosis`),  # If Date of Diagnosis is valid
        as.integer(round(as.numeric(difftime(as.Date(`Date of Diagnosis`, "%m/%d/%Y"), as.Date(DOB, "%m/%d/%Y"), units = "days")) / 365.25)),
        NA_integer_
      )
    )
  )

# 4. Clean `Final Pathology` column 
mst_list <- mst_list %>%
  mutate(`Final Pathology` = clean_pathology(`Final Pathology`))

# 5. Clean `Surgery` column 
mst_list <- mst_list %>%
  mutate(SURGERY = clean_surgery(SURGERY))


# 6. Clean `Res_SEQ-WTS (Y/N)`
mst_list <- mst_list %>%
  mutate(`Res_SEQ-WTS (Y/N)` = ifelse(`Res_SEQ-WTS (Y/N)` == "Y", "Y", "N"))

# 7. Drop DOB
mst_list <- mst_list %>% dplyr::select(-DOB)

# 8. Drop columns not needed for analysis 
mst_list <- mst_list %>% dplyr::select(-`Res_DNA (CARIS) (Y/N)`, -`Res_RNA (CARIS) (Y/N)`) # nolint: line_length_linter.