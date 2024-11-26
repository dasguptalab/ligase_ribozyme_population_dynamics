#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# round numbers
rounds <- c(1, 2, 3, 4, 5, 6, 7, 8)

# % diversity per round
#diversity_doped <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20, 97.30, 92.40, 86.43)
diversity <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20)

# read in cluster identity data
r8_peaks_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_cluster_identity/07a_clustered/r8_S8_L001_formatted_above9_cluster_peaks_identity_table.csv")
peaks_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_cluster_identity/07a_clustered/peaks_identity_table.csv")
#seqs_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_cluster_identity/07a_clustered/sequences_identity_table.csv")

# TO-DO: double check for duplicate data... (low ID cluster peak entries)
# read in cluster family sequence data
r8_seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_cluster_families/07a_clustered/r8_S8_L001_cluster_families_table.csv")
#r8_peaks_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_cluster_families/07a_clustered/r8_cluster_peaks_families_table.csv")
r8_peaks_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_cluster_families/07a_clustered/r8_S8_L001_cluster_peaks_families_table.csv")

# setup counts data frame
#seq_log_counts <- data.frame(
#  run_name = NA,
#  seq_num = NA,
#  log_counts = NA
#)

# loop over each run
#for (run_num in 1:8) {
  # TO-DO: change to combined seq data files
  # setup the file name
#  file_in <- paste("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_cluster_families/07a_clustered/r", run_name, "_cluster_peaks_table.csv", sep="")
  # read in the counts for the top 10 sequences of each round, per round
#  seqs_family_top10 <- read.csv(file_in, nrows=10)
  # loop over the top 10 seqs
#  for (seq_num in 1:10) {
    # set tmp index
#    tmp_index <- run_num-1
    # set the index
#    index <- run_num+(tmp_index*8)
    # setup the column name
#    curr_col <- paste("r", run_num, "_counts", sep="")
    # calculate the log counts for the top 10 sequences of each round, per round
    
#    }
#}

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# list of cluster IDs
cluster_list <- unique(r8_seqs_family$cluster_ID)
cluster_list_out <- cluster_list+1

# setup data frame length
data_length <- (max(cluster_list)+1)*8

# data frame of cluster abundances
cluster_abundances <- data.frame(
  cluster_ID = rep(NA, data_length), 
  run_name = rep(NA, data_length),
  cluster_color = rep(NA, data_length),
  frac_abundance = rep(NA, data_length)
)

# calculate fraction abundance per round
for (cluster_num in 0:max(cluster_list)) {
  # loop over each run
  for (run_num in 1:8) {
    # set the index
    index <- run_num+(cluster_num*8)
    # update family number for publishing
    cluster_out <- cluster_num+1
    # add cluster ID
    cluster_abundances$cluster_ID[index] <- paste("family", cluster_out, sep="_")
    # add run name
    cluster_abundances$run_name[index] <- run_num
    # add cluster plotting color
    cluster_abundances$cluster_color[index] <- safe_colors[cluster_out]
    # setup the column name
    curr_col <- paste("r", run_num, "_counts", sep="")
    # add fraction abundance
    cluster_abundances$frac_abundance[index] <- sum(r8_seqs_family[r8_seqs_family$cluster_ID == cluster_num, curr_col])/quality[1]
  }
}

# calculate % unique per round
unique_reads <- diversity

# create data frame for plotting
unique_rounds <- data.frame(
  round_nums = rounds,
  unique_percent = unique_reads
)

# line plot with % of unique sequences per round
ggplot(data=unique_rounds, aes(round_nums, unique_percent)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylab("Percent Unique Sequences") +
  xlab("Round Number")

# line plots with round 8 clustering identity info
ggplot(data=r8_peaks_identity, aes(x = cluster_ID)) +
  geom_line(aes(y = avg_ID, color = "Avg")) +
  #geom_line(aes(y = highest_ID, color = "Highest")) +
  geom_line(aes(y = lowest_ID, color = "Lowest")) +
  geom_point(aes(y = avg_ID, color = "Avg")) +
  #geom_point(aes(y = highest_ID, color = "Highest")) +
  geom_point(aes(y = lowest_ID, color = "Lowest")) +
  scale_color_manual(name = "Statistic", values = c("Avg" = safe_colors[5], "Highest" = safe_colors[9], "Lowest" = safe_colors[4])) +
  theme_bw() +
  scale_x_continuous(breaks = breaks_pretty()) +
  ylab("Percent Identity to Peak") +
  xlab("Cluster Number")

# line plot with fraction abundance per round for each of the top 10 families
ggplot(data=cluster_abundances, aes(x=run_name, y=frac_abundance, group=cluster_ID, color=cluster_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(labels = cluster_list_out, breaks = safe_colors, name = "Family", guide = "legend") +
  ylab("Fraction Abundance") +
  xlab("Round Number")

# round 8 peak seqs
# setup length of counts
counts_length <- 13*8
# tmp data frame for peak counts
subset_peak_counts <- data.frame(
  round_num = rep(NA, counts_length),
  seq_num = rep(NA, counts_length),
  counts = rep(NA, counts_length)
)
# data frame for peak counts
peak_counts <- data.frame()
# seq numbers for plotting
for (round_index in 1:8) {
  # set the index
  index <- ((round_index-1)*8)+1
  # setup the column name
  curr_col <- paste("r", round_index, "_counts", sep="")
  # add round num
  subset_peak_counts$round_num <- rep(round_index, 13)
  # add seq num
  subset_peak_counts$seq_num <- seq(from = 1, to = 13, by = 1)
  # add peak seq counts
  subset_peak_counts$counts <- r8_peaks_family[, curr_col]
  # add rows
  peak_counts <- rbind(peak_counts, subset_peak_counts)
}

# heatmaps with the log counts for each of the top peak sequences from round 8
ggplot(data = peak_counts, aes(round_num, seq_num, fill= log(counts))) + 
  theme_bw() +
  geom_tile() +
  ylab("Sequence Number") +
  xlab("Round Number")

# heatmaps with the log counts for each of the top 10 sequences per round


# heatmaps with the log counts for each of all sequences per round

