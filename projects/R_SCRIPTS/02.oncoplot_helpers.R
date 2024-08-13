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
            show_column_names = FALSE,
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

git clone https://<USERNAME>:<PASSWORD>@github.com/<USERNAME>/<REPO_NAME>.git
git clone "https://kasondechw:zxhL7%0f`T*3D<i78y8`L9g8S$D!@github.com/mxr893/dang_lab.git"