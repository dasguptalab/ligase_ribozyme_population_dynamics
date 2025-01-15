#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# install logo tool
# https://omarwagih.github.io/ggseqlogo/
#devtools::install_github("omarwagih/ggseqlogo")

# import libraries
library(ggplot2)
#library(scales)
library(rcartocolor)
library(stringr)
require(ggseqlogo)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/04a_family_sequence_logos"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# read in cluster sequence family data
seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/r8_S8_L001_formatted_above9_cluster_sequences_identity_table_atLeast90.csv")

# convert Ts to Us
seqs_input$sequence <- gsub("T", "U", seqs_input$sequence)

# read in sequences that have at least 90% identity to any peak
seqs_identities <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/00a_family_identification_above9/family_identities_above9_atLeast90.csv")

# subset to keep round 8 data
seqs_identities <- seqs_identities[seqs_identities$run_name == 8,]

# initialize data frame
seqs_family <- data.frame()
peak_cluster_IDs <- NULL

# loop over each sequences that has at least 90% identity to any peak
for (seq_num in 1:nrow(seqs_identities)) {
  # keep count data for sequences that have at least 90% identity to any peak
  seqs_90_data <- seqs_input[seqs_input$sequence_ID == strsplit(seqs_identities$sequence_ID[seq_num], "_")[[1]][2],]
  # add counts
  seqs_family <- rbind(seqs_family, seqs_90_data)
  # add peak cluster ID to vector
  peak_cluster_IDs <- c(peak_cluster_IDs, rep(seqs_identities$peak_cluster_ID[seq_num], nrow(seqs_90_data)))
}

# add peak cluster IDs
seqs_family <- cbind(seqs_family, peak_cluster_IDs)

# list of cluster IDs
fam_list_out <- seq(1, 13)

# list of cluster IDs in order of abundance in round 8
r8_fams <- data.frame(
  fam_ID = fam_list_out,
  cluster_ID = c(1, 3, 0, 2, 5, 8, 4, 7, 11, 6, 10, 12, 9)
)

# list of cluster IDs
cluster_list <- unique(seqs_family$cluster_ID)

# loop over each cluster
for (cluster_num in min(cluster_list):max(cluster_list)) {
  #cluster_num <- 1
  # convert list of sequences into a matrix
  seqs_matrix <- do.call(rbind, type.convert(strsplit(seqs_family[seqs_family$cluster_ID == cluster_num, "sequence"], ""), as.is = TRUE))
  # set family number
  fam_num <- r8_fams[r8_fams$cluster_ID == cluster_num, "fam_ID"]
  # set family plot title
  run_title <- paste("Family", fam_num, "Sequence Logo")
  # create sequence logo
  logo_plot <- ggplot() + 
    geom_logo(seqs_family[seqs_family$cluster_ID == cluster_num, "sequence"], method = 'prob') + 
    # left stem
    annotate('text', x=9, y=1.05, label='Left Stem') + 
    annotate('rect', xmin = 5.5, xmax = 12.5, ymin = -0.01, ymax = 1.01, 
             alpha = .1, col=safe_colors[2], fill=safe_colors[3], linewidth = 1) + 
    # right stem
    annotate('text', x=29, y=1.05, label='Right Stem') + 
    annotate('rect', xmin = 25.5, xmax = 32.5, ymin = -0.01, ymax = 1.01, 
             alpha = .1, col=safe_colors[2], fill=safe_colors[3], linewidth = 1) + 
    # overhang
    annotate('text', x=18, y=1.05, label='Overhang') + 
    annotate('rect', xmin = 15.5, xmax = 20.5, ymin = -0.01, ymax = 1.01, 
             alpha = .1, col=safe_colors[1], fill=safe_colors[3], linewidth = 1) + 
    theme_logo() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
    ggtitle(run_title) +
    theme(plot.title = element_text(hjust = 0.5))
  # save the plot
  exportFile <- paste(out_dir, "/family", fam_num, "_sequence_logo.png", sep = "")
  png(exportFile, units="in", width=10, height=4, res=300)
  print(logo_plot)
  dev.off()
}

# export plotting data
write.csv(seqs_family, file = paste(out_dir, "/family_sequences.csv", sep = ""), row.names = FALSE, quote = FALSE)
