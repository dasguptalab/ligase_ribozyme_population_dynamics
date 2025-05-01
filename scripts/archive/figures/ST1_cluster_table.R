#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# import libraries
library(ggplot2)
library(scales)
library(rcartocolor)
#library(plyr)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/ST1_cluster_table"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
unique <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)

# read in cluster family sequence data
#seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized/07a_clustered/cluster_peaks_table.csv")
#r8_seqs_family <- seqs_family[seqs_family$run_name == "r8",]
r8_seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv")

# set the round num
roundNum <- 8

# list of cluster IDs
#cluster_list <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
cluster_list <- r8_seqs_family$cluster_ID

# color blind safe plotting palette
#safe_colors <- c(carto_pal(name="Safe"), "#000000")
safe_colors <- carto_pal(name="Safe")[1:length(cluster_list)]

# create cluster data frame
cluster_data <- data.frame(
  cluster_ID = cluster_list,
  cluster_color = safe_colors
)

# add peak sequences
cluster_data_out <- merge(cluster_data, r8_seqs_family, by = "cluster_ID")

# add abundances
cluster_data_out$read_abun <- 100*cluster_data_out$read_counts/quality[roundNum]

# sort cluster data
cluster_data_out <- cluster_data_out[order(cluster_data_out$read_abun, decreasing = TRUE),]  

# add family numbers
#cluster_data_out$fam_num <- seq(from = 1, to = length(cluster_list), by = 1)
  
# export data
write.csv(cluster_data_out, file = paste(out_dir, "/r8_cluster_count_data.csv", sep = ""), row.names = FALSE, quote = FALSE)

