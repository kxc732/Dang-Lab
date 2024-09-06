#!/bin/R Scripts 
# Date: 8/13/2024  
# File: 02.oncoplot_helpers.R
# Author: Kasonde Chewe 
# > Description: Helper functions use to pre-process data into correct format and generate summary oncoplot along with additional 
# >              mutation summary report for cohort size N in JMP PAC datasets

# Required Libraries 
library(ComplexHeatmap)


# Function to generate oncoplot, adjust font sizes and visual params 
# See Additional Params Below 
generate_oncoplot <- function(matrix, title = NULL) {
  oncoPrint(matrix,
            alter_fun = alter_fun_list,
            alter_fun_is_vectorized = FALSE,
            col = consequence_colors,
            column_title = title,
            column_title_gp = gpar(fontsize = 24),  # Adjust column title font size
            show_column_names = TRUE,
            show_row_names = TRUE,
            row_names_gp = gpar(fontsize = 16),  # Adjust row names font size
            heatmap_legend_param = list(
              title = "Alterations", 
              title_gp = gpar(fontsize = 16),  # Adjust legend title font size
              at = names(consequence_colors), 
              labels = names(consequence_colors)
            ),
            row_order = NULL,
            remove_empty_rows = FALSE,
            show_pct = TRUE,  # Show mutation frequencies (percentages)
            pct_gp = gpar(fontsize = 16)  # Adjust font size for percentages
  )
}


# Add missing mutation types to the colors vector
consequence_colors <- c(
  "FRAMESHIFT" = "green", 
  "MISSENSE" = "blue", 
  "SPLICING" = "orange", 
  "NONSENSE" = "purple", 
  "CODON INSERTION" = "yellow", 
  "START LOST" = "brown", 
  "GAIN" = "red", 
  "LOSS" = "black",
  "CODON DELETION" = "pink",  # Add color for CODON DELETION
  "SYNONYMOUS" = "lightblue", # Add color for SYNONYMOUS
  "NONCODING" = "#008080",       # Add color for NONCODING
  "CODON_CHANGE_PLUS_CODON_DELETION" = "magenta", # Add color for CODON_CHANGE_PLUS_CODON_DELETION
  "SILENT" = "cyan",          # Add color for SILENT
  "UTR" = "lightgreen",       # Add color for UTR
  "NA" = "white",             # Ensure color for NA
  "Normal" = "#CCCCCC"        # Ensure color for Normal
)

# Update the alter_fun_list to include all necessary functions
alter_fun_list <- list(
  Normal = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["NORMAL"], col = NA))
  },
  FRAMESHIFT = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["FRAMESHIFT"], col = NA))
  },
  MISSENSE = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["MISSENSE"], col = NA))
  },
  SPLICING = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["SPLICING"], col = NA))
  },
  NONSENSE = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["NONSENSE"], col = NA))
  },
  `CODON INSERTION` = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["CODON INSERTION"], col = NA))
  },
  `START LOST` = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["START LOST"], col = NA))
  },
  GAIN = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = NA, col = "red", lwd = 3))
  },
  LOSS = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = NA, col = "black", lwd = 3))
  },
  `CODON DELETION` = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["CODON DELETION"], col = NA))
  },
  SYNONYMOUS = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["SYNONYMOUS"], col = NA))
  },
  NONCODING = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["NONCODING"], col = NA))
  },
  `CODON CHANGE PLUS CODON DELETION` = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["CODON CHANGE PLUS CODON DELETION"], col = NA))
  },
  SILENT = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["SILENT"], col = NA))
  },
  UTR = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["UTR"], col = NA))
  },
  "NA" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = consequence_colors["NA"], col = NA))
  },
  background = function(x, y, w, h){
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA))
  }
)

# Define a function to categorize obesity based on BMI
categorize_obesity <- function(bmi) {
  ifelse(as.numeric(bmi) >= 30, "Y", "N")
}


# # Update colors 
# # Update ne
# # Create the bottom annotation using the risk_factors dataframe
# bottom_annotation <- HeatmapAnnotation(
#   df = risk_factors[, -1], # Exclude the JMP_PAC ID
#   col = list(
#     Diabetes = c("Y" = "red", "N" = "blue"),
#     Pancreatitis1 = c("Y" = "orange", "N" = "green"),
#     Pancreatitis2 = c("chronic" = "purple", "acute" = "pink", "NA" = "gray", "unknown"="black"),
#     Hyperlipidemia = c("Y" = "yellow", "N" = "black"),
#     Fatty_Liver_Disease = c("Y" = "cyan", "N" = "brown", "NA" = "lightgray"), # Add color for NA
#     Living_Status = c("Y" = "#32CD32", "N" = "navy"),
#     Hypertension = c("Y" = "gold", "N" = "#C0C0C0"),
#     Obesity = c("Y" = "darkred", "N" = "darkblue")
#   ),
#   annotation_height = unit(1, "cm"),
#   border = TRUE # Add a border to the annotations
# )

# Generate oncoplot with bottom annotation
generate_oncoplot_with_anno <- function(matrix, title = NULL) {
  oncoPrint(matrix,
            alter_fun = alter_fun_list,
            alter_fun_is_vectorized = FALSE,
            col = consequence_colors,
            column_title = title,
            column_title_gp = gpar(fontsize = 24),
            show_column_names = TRUE,
            show_row_names = TRUE,
            row_names_gp = gpar(fontsize = 16),
            heatmap_legend_param = list(
              title = "Alterations", 
              title_gp = gpar(fontsize = 16),
              at = names(consequence_colors), 
              labels = names(consequence_colors)
            ),
            row_order = NULL,
            remove_empty_rows = FALSE,
            show_pct = TRUE,
            pct_gp = gpar(fontsize = 16),
            bottom_annotation = bottom_annotation # Add the bottom annotation here
  )
}

# in R code 
# function must be moved to other helper section 
# function accepts a dataframe as param
# second param is will be the column name "ID"
# third param is the filepath column name "path"
# the function must read reach tsv file and create a gene counts matrix for all 120 patients 
# the function must then return this matrix 

# Load necessary libraries
library(dplyr)

# Define the helper function
create_gene_counts_matrix <- function(df, id_col, path_col) {
  
  # Initialize an empty list to store the gene count data frames
  gene_counts_list <- list()
  
  # Iterate over each row of the dataframe
  for (i in 1:nrow(df)) {
    # Extract the ID and filepath for the current patient
    patient_id <- df[[id_col]][i]
    file_path <- as.character(df[[path_col]][i])
    
    # Read the TSV file into a data frame, without header and using default column names
    gene_counts <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    # Ensure the first column is the gene names
    gene_names <- gene_counts[, 1]
    
    # The second column contains the counts
    counts <- gene_counts[, 2]
    
    # Create a named vector for counts
    named_counts <- setNames(counts, gene_names)
    
    # Store the named vector in the list, using the patient ID as the list name
    gene_counts_list[[patient_id]] <- named_counts
  }
  
  # Combine all the named vectors into a data frame, filling missing values with NA
  gene_counts_matrix <- bind_rows(gene_counts_list, .id = "ID") %>%
                        column_to_rownames(var = "ID") %>%
                        t()  # Transpose the matrix so IDs are columns and genes are rows
  
  # Return the resulting gene counts matrix
  return(gene_counts_matrix)
}

# Define the helper function to extract gene values
extract_gene_values <- function(gene_counts_matrix, gene_name) {
  
  # Check if the gene exists in the matrix
  if (!gene_name %in% rownames(gene_counts_matrix)) {
    stop(paste("Gene", gene_name, "not found in the matrix."))
  }
  
  # Extract the values for the specified gene across all columns (patients)
  gene_values <- gene_counts_matrix[gene_name, ]
  
  # Convert the gene values to a data frame
  gene_exp_df <- as.data.frame(gene_values)
  
  # Add patient IDs as a new column
  gene_exp_df$`JMP_PAC ID`<- rownames(gene_exp_df)
  
  # Reset the row names
  rownames(gene_exp_df) <- NULL
  
  # Rename the gene expression column to 'Expression (TPM)'
  colnames(gene_exp_df)[colnames(gene_exp_df) == "gene_values"] <- "Expression (TPM)"
  
  # Move the ID column to the first position
  gene_exp_df <- gene_exp_df %>% dplyr::select(`JMP_PAC ID`, everything())
  
  # Return the final data frame
  return(gene_exp_df)
}

