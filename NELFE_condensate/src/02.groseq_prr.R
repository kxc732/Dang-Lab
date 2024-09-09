# Title: PRR Analysis for GRO-seq Data
# Author(s): 
# Associated Manuscript: [Title of Manuscript]
# Journal: [Journal Name], 2024
# Corresponding Author: 
# Date: [YYYY-MM-DD]
# Purpose: This script computes GRO-seq Pausing indices (PI) and promoter release ratios (PRR) 
# for gene expression analysis using GRO-seq data.

# Libraries 
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(BRGenomics)
library(genomation)
library(ggplot2)
library(org.Hs.eg.db)
library(biomaRt)
library(Homo.sapiens)
library(ggpointdensity)
library(AnnotationDbi)
library(ggsignif)
library(dplyr)
library(writexl)

# Helper functions in groseq_helpers.R
source("helpers/groseq_helpers.R")

# GRO-seq Coverage 
ev_bw_path <- "gro-ev.p.igdup.bw"
ko_bw_path <- "gro-ko.p.igdup.bw"

ev_bw <- keepSeqlevels(import(ev_bw_path), paste0("chr", c(1:22)), pruning.mode = "coarse")
ko_bw <- keepSeqlevels(import(ko_bw_path), paste0("chr", c(1:22)), pruning.mode = "coarse")

# Transcript Annotation
TxDb(Homo.sapiens) <- keepSeqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene, paste0("chr", c(1:22)), pruning.mode = "coarse")
tx <- transcriptsBy(Homo.sapiens, columns = "SYMBOL")
gr <- keepSeqlevels(unlist(tx), paste0("chr", c(1:22)), pruning.mode = "coarse")

# Overlap ATAC, H3K27ac ChIP-seq and Gro-seq
ah_11356 <- read.delim("ah.11536.bed", sep = "\t", header = FALSE, col.names = c("seqnames", "start", "end"))
ah_11356_gr <- GRanges(seqnames = ah_11356$seqnames, ranges = IRanges(start = ah_11356$start + 1, end = ah_11356$end))
overlaps <- findOverlaps(gr, ah_11356_gr)

overlapping_gr <- gr[queryHits(overlaps)]
gr <- overlapping_gr

# Clean up
rm(ah_11356, ah_11356_gr, overlaps)

# Step 4: Compute Gene Density in Gene Body and Promoter Regions
pr <- genebodies(gr, start = -250, end = 300, fix.end = "start")
gb <- genebodies(gr, start = 301, end = 2000, fix.end = "start")

# Step 5: Compute Pausing Index (PI) Using BRGenomics
pidx_ev <- getPausingIndices(ev_bw, pr, gb)
pidx_ko <- getPausingIndices(ko_bw, pr, gb)

# Step 6: Create DataFrame for PI Values
pausing_data <- data.frame(
  gene = pr$SYMBOL,
  PI_EV = pidx_ev,
  PI_KO = pidx_ko
)

# Clear env 
rm(pidx_ev, pidx_ko, ev_bw, ko_bw, pr, gb)

# Step 7: Handle Infinite and Zero PI Values
pausing_data$PI_EV <- replace(pausing_data$PI_EV, is.infinite(pausing_data$PI_EV), NA)
pausing_data$PI_KO <- replace(pausing_data$PI_KO, is.infinite(pausing_data$PI_KO), NA)

# Ensure data is numeric
pausing_data$PI_EV <- as.numeric(pausing_data$PI_EV)
pausing_data$PI_KO <- as.numeric(pausing_data$PI_KO)

# Step 8: Aggregate and De-duplicate
pausing_data_median <- aggregate(cbind(PI_EV, PI_KO) ~ gene, data = pausing_data, FUN = median, na.rm = TRUE)
pausing_data_median$actively_paused_PI_EV <- pausing_data_median$PI_EV > 2
pausing_data_median$actively_paused_PI_KO <- pausing_data_median$PI_KO > 2

# Step 9: Compute Promoter Release Ratio (PRR)
pausing_data_median$PRR_EV <- 1 / pausing_data_median$PI_EV
pausing_data_median$PRR_KO <- 1 / pausing_data_median$PI_KO

# Step 10: IQR Filtering

pausing_data_median <- iqr_filter(pausing_data_median, 'PI_EV', 'PI_KO')
pausing_data_median <- iqr_filter(pausing_data_median, 'PRR_EV', 'PRR_KO')

# Log2 transformation for plotting
copy <- pausing_data_median
copy$PI_EV <- log2(copy$PI_EV) + 0.0001
copy$PI_KO <- log2(copy$PI_KO) + 0.0001
copy$PRR_EV <- log2(copy$PRR_EV) + 0.0001
copy$PRR_KO <- log2(copy$PRR_KO) + 0.0001

# Step 11: Generate Boxplots for PI and PRR
boxplot_c2(data = copy, col1 = 'PI_EV', col2 = 'PI_KO')
boxplot_c2(data = copy, col1 = 'PRR_EV', col2 = 'PRR_KO')

# Clean up
rm(copy)

# Step 12: Compute log2FC for PI and PRR
pausing_data_median$log2FC_PI <- log2(pausing_data_median$PI_KO + 0.0001) - log2(pausing_data_median$PI_EV + 0.0001)
pausing_data_median$log2FC_PRR <- log2(pausing_data_median$PRR_KO + 0.0001) - log2(pausing_data_median$PRR_EV + 0.0001)

# Step 13: Summarize Data Using Tertiles
pausing_data_median$tertile_PRR <- cut(
  pausing_data_median$log2FC_PRR,
  breaks = c(-Inf, 1, 4, Inf),
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)

# Step 14: Save Data to Excel Files
writexl::write_xlsx(pausing_data_median, "PI_summary.xlsx")

# Optional: Perform pairwise Wilcoxon test on tertile data and generate boxplots for log2FC PRR
pairwise_results <- pairwise.wilcox.test(pausing_data_median$log2FC_PRR, pausing_data_median$tertile_PRR, p.adjust.method = "bonferroni")

pdf("tertile_boxplots.pdf", width = 10, height = 10)
ggplot(pausing_data_median, aes(x = tertile_PRR, y = log2FC_PRR, fill = tertile_PRR, color = tertile_PRR)) +
  geom_boxplot(position = position_dodge(0.8), fill = NA, lwd = 1.0, outlier.shape = NA) +
  theme_bw() +
  geom_signif(comparisons = combn(levels(pausing_data_median$tertile_PRR), 2, simplify = FALSE),
              map_signif_level = TRUE,
              textsize = 4,
              y_position = c(8, 9, 10)) +  
  scale_color_manual(values = c("blue", "black", "red")) +
  labs(title = "Tertile Boxplot of log2FC PRR (KO/EV)",
       x = "Tertile PRR",
       y = "log2 Fold-Change PRR") +
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"))
dev.off()

# Print summary statistics
tertile_counts <- table(pausing_data_median$tertile_PRR)
highly_elongated_genes <- tertile_counts["High"]
moderately_elongated_genes <- tertile_counts["Medium"]
not_actively_elongated_genes <- tertile_counts["Low"]

cat("Highly elongated genes =", highly_elongated_genes, "\n")
cat("Moderately elongated genes =", moderately_elongated_genes, "\n")
cat("Not actively elongated genes =", not_actively_elongated_genes, "\n")


