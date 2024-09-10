#!/usr/bin/env Rscript
# Title: PRR Analysis for GRO-seq Data Helper Functions
# Author(s): 
# Associated Manuscript: [Title of Manuscript]
# Journal: [Journal Name], 2024
# Corresponding Author: 
# Date: [YYYY-MM-DD]
# Purpose: This script includes associated helper functions in GRO-seq PRR and PI 


# Function for creating a density correlatin plot for two columns
correlation_c2 <- function(data, col1, col2, title=NULL){
  # Check if the specified columns exist in the dataframe
  if (!col1 %in% names(data) | !col2 %in% names(data)) {
    stop("Specified columns do not exist in the dataframe")
  }
}

# Function for creating boxplots for two columns
boxplot_c2 <- function(data, col1, col2, title=NULL) {
  # Check if the specified columns exist in the dataframe
  if (!col1 %in% names(data) | !col2 %in% names(data)) {
    stop("Specified columns do not exist in the dataframe")
  }
  
  # Reshape the data from wide to long format for ggplot
  data_long <- reshape2::melt(data, measure.vars = c(col1, col2))
  
  # Create the box plot
  p <- ggplot(data_long, aes(x = variable, y = value, fill= variable, color = variable)) + 
    theme_bw() +
    geom_boxplot(outlier.color = "red", outlier.shape = NA, fill = NA) +
    scale_fill_manual(values = c("red", "black")) + 
    scale_x_discrete(labels = c(col1, col2)) +
    labs(title = paste("Box plot for", col1, "and", col2, title),
         x = "Condition",
         y = "Value") 
  # Print the plot
  print(p)
}


# Function for the IQR filtering function
iqr_filter <- function(data, col1, col2) {
  # Calculate the IQR for col1
  Q1_col1 <- quantile(data[[col1]], 0.05, na.rm = TRUE)
  Q3_col1 <- quantile(data[[col1]], 0.95, na.rm = TRUE)
  IQR_col1 <- Q3_col1 - Q1_col1
  
  # Calculate the IQR for col2
  Q1_col2 <- quantile(data[[col2]], 0.05, na.rm = TRUE)
  Q3_col2 <- quantile(data[[col2]], 0.95, na.rm = TRUE)
  IQR_col2 <- Q3_col2 - Q1_col2
  
  # Function for the cutoff for outliers
  cutoff_low_col1 <- Q1_col1 - 1.5 * IQR_col1
  cutoff_high_col1 <- Q3_col1 + 1.5 * IQR_col1
  cutoff_low_col2 <- Q1_col2 - 1.5 * IQR_col2
  cutoff_high_col2 <- Q3_col2 + 1.5 * IQR_col2
  
  # Filter the data based on the IQR cutoff
  data_filtered <- data[data[[col1]] > cutoff_low_col1 & data[[col1]] < cutoff_high_col1 & 
                          data[[col2]] > cutoff_low_col2 & data[[col2]] < cutoff_high_col2, ]
  
  return(data_filtered)
}


# Function to generate ECDF from 2 column data 
generate_ecdf <- function(df, col1, col2, ps = 1){
  # Compute log2
  df$log2_EV <- log2(df$col1 + ps)
  df$log2_KO <- log2(df$col2 + ps)
  ecdf_ev <- ecdf(df$log2_EV)
  ecdf_ko <- ecdf(df$log2_KO)
  
  # Plot the ECDFs for log2-transformed PI_EV and PI_KO with pseudocount
  plot(ecdf_ev, col = "blue", main = "ECDF of log2(PI_EV) and log2(PI_KO)", 
       xlab = "log2(PI + pseudocount) Value", ylab = "ECDF", xlim = c(-8, 6))
  lines(ecdf_ko, col = "red")
  
  # Add a legend
  legend("bottomright", legend = c("log2(PI_EV)", "log2(PI_KO)"), col = c("blue", "red"), lty = 1)
}

# Function creates 3 low, medium and high tertiled boxplots on a single graph
create_double_tertile_boxplot<- function(data, ev_col_name, ko_col_name, title=NULL, saveas=NULL, csv=NULL) {
  ev_col_name <- as.character(substitute(ev_col_name))
  ko_col_name <- as.character(substitute(ko_col_name))
  
  # Add a column for log fold change (log2(KO/EV))
  data$logFC <- log2(data[[ko_col_name]] / data[[ev_col_name]] + 1)  # Adding 1 to prevent division by zero
  
  
  
  # Remove rows with NA in logFC
  data <- na.omit(data, cols = "logFC")
  
  # Define criteria for high, medium, and low tertiles
  data$Tertile <- cut(data$logFC,
                      breaks = c(-Inf, 1.4, 3, Inf),
                      labels = c("Low", "Medium", "High"),
                      include.lowest = TRUE)
  
  View(data)
  
  # save data to a csv 
  if (!is.null(csv)) {
    write.csv(data, csv, row.names = FALSE)}
  
  
  # Create the boxplot for log fold change
  p <- ggplot(data, aes(x = Tertile, y = logFC, fill = Tertile, color = Tertile)) +
    theme_bw() +
    geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, fill=NA) +
    labs(x = 'Tertile', y = 'Log2FC(KO/EV)', title = title) +
    scale_color_manual(values = c("blue", "green", "red")) +
    ylim(0, 6.5) +
    theme(axis.text.x  = element_text(size = 20),   # Adjust x-axis label font size
          axis.text.y = element_text(size = 20),    # Adjust y-axis label font size
          axis.title.x  = element_text(size = 20),  # Adjust x-axis title font size
          axis.title.y = element_text(size = 20))   # Adjust y-axis title font size
  
  # Save or print the plot
  if (!is.null(saveas)) {
    ggsave(saveas, p, width = 10, height = 10, dpi = 600, device = "pdf")
  } else {
    print(p)
  }
} 


# Function to compute KO, EV log fold change
generate_pausing_ratio_across_tertiles <- function(data, ev_col_name, ko_col_name, pseudocount = 1, title=NULL, saveas=NULL) {
  if (any(is.na(data[[ev_col_name]])) | any(is.na(data[[ko_col_name]]))) {
    stop("Data contains missing values in critical columns. Please handle them before proceeding.")
  }
  
  # Offset the data to ensure all values are positive
  min_value <- min(c(min(data[[ev_col_name]], na.rm = TRUE), min(data[[ko_col_name]], na.rm = TRUE)))
  offset <- abs(min_value) + 1
  data[[ev_col_name]] <- data[[ev_col_name]] + offset
  data[[ko_col_name]] <- data[[ko_col_name]] + offset
  
  # Log transformation after offsetting
  data[[ev_col_name]] <- log2(data[[ev_col_name]] + pseudocount)
  data[[ko_col_name]] <- log2(data[[ko_col_name]] + pseudocount)
  
  # Calculate log fold change
  data$logFC <- data[[ko_col_name]] - data[[ev_col_name]]
  
  
  # Define tertiles for log fold change
  data$Tertile <- cut(data$logFC,
                      breaks = quantile(data$logFC, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      labels = c("Low", "Medium", "High"),
                      include.lowest = TRUE)
  
  # Convert Tertile to a factor with levels in the order you want
  data$Tertile <- factor(data$Tertile, levels = c("Low", "Medium", "High"), ordered = TRUE)
  
  # Extract the Transcript IDs for each tertile
  data$TranscriptID <- data[,2]
  
  # Filter the original data for each tertile category
  low_data <- data[data$Tertile == "Low", c(ev_col_name, ko_col_name)]
  medium_data <- data[data$Tertile == "Medium", c(ev_col_name, ko_col_name)]
  high_data <- data[data$Tertile == "High", c(ev_col_name, ko_col_name)]
  
  # Combine filtered data into one dataframe with a new 'Tertile' column for plotting
  combined_data <- rbind(
    cbind(low_data, Tertile = 'Low'),
    cbind(medium_data, Tertile = 'Medium'),
    cbind(high_data, Tertile = 'High')
  )
  # Reshape the combined data from wide to long format for plotting
  long_combined_data <- combined_data %>%
    pivot_longer(cols = c(ev_col_name, ko_col_name), names_to = "Condition", values_to = "Value") %>%
    mutate(Tertile = factor(Tertile, levels = c("Low", "Medium", "High"), ordered = TRUE))
  
  # Create the boxplot for original EV and KO values for each tertile
  p <- ggplot(long_combined_data, aes(x = Tertile, y = Value, fill = Condition)) +
    geom_boxplot(aes(color = Condition), position = position_dodge(0.8), outlier.shape = NA) +
    theme_minimal() +
    labs(x = 'Tertile', y = 'Value', title = title) +
    scale_fill_manual(values = c("white", "white")) +
    scale_color_manual(values = c("red", "black")) +
    theme(legend.position = "bottom")
  
  # Save or print the plot
  if (!is.null(saveas)) {
    ggsave(saveas, p, width = 12, height = 10, dpi = 600, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
  } else {
    print(p)
  }
  
  # Return a list of genes in each tertile
  return(list(Low = data$TranscriptID[data$Tertile == "Low"],
              Medium = data$TranscriptID[data$Tertile == "Medium"],
              High = data$TranscriptID[data$Tertile == "High"]))
} 


add_genomic_coordinates <- function(df, genes, all_regions) {
  # Initialize connection to Ensembl
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Retrieve genomic coordinates
  gene_coords <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position', 'strand'),
                       filters = 'hgnc_symbol',
                       values = genes,
                       mart = ensembl)
  
  # Handle multiple entries per gene by selecting the first entry per gene symbol
  gene_coords <- gene_coords[!duplicated(gene_coords$hgnc_symbol), ]
  
  # Rename columns for BED-like format
  names(gene_coords) <- c("gene.value", "chr", "start", "end", "strand")
  
  # Edit chromosome values to include 'chr' prefix if numeric
  gene_coords$chr <- ifelse(grepl("^[0-9]+$", gene_coords$chr), paste0("chr", gene_coords$chr), gene_coords$chr)
  
  # Normalize strand values to '+' or '-'
  gene_coords$strand <- ifelse(gene_coords$strand == 1, "+", ifelse(gene_coords$strand == -1, "-", gene_coords$strand))
  
  # Merge with the original dataframe
  df <- merge(df, gene_coords, by = "gene.value", all.x = TRUE)
  
  # Identify entries with non-standard chromosome names
  non_standard_genes <- df %>% filter(!grepl("^chr[0-9XY]+$", chr))
  
  # Correct non-standard entries using all_regions data
  for(i in seq_len(nrow(non_standard_genes))) {
    gene_name <- non_standard_genes$gene.value[i]
    correct_chr <- all_regions %>% filter(Gene.Name == gene_name) %>% pull(Chr) %>% unique()
    if(length(correct_chr) == 1) {
      df$chr[df$gene.value == gene_name] <- correct_chr
    }
  }
  
  # Standardize strand notation in the dataframe
  df$strand <- ifelse(df$strand == 1, "+", ifelse(df$strand == -1, "-", df$strand))
  
  # Drop rows where chr is NA
  df <- df[!is.na(df$chr), ]
  
  # Drop rows where start is NA
  df <- df[!is.na(df$start), ]
  
  # Reorder the DataFrame so that genomic columns are the first four after 'gene.value'
  df <- df %>%
    dplyr::select(gene.value, chr, start, end, strand, everything())
  
  # Return the updated dataframe
  return(df)
}
