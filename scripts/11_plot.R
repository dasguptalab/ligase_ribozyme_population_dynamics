#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)

# color blind safe plotting palette
safe_colors <- carto_pal(name="Safe")

# round numbers
rounds <- c(1, 2, 3, 4, 5, 6, 7, 8)

# % diversity per round
#diversity_doped <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20, 97.30, 92.40, 86.43)
diversity <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20)

# read in cluster identity data
peaks_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_cluster_identity/07a_clustered/cluster_peaks_identity_table.csv")
seqs_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_cluster_identity/07a_clustered/cluster_sequences_identity_table.csv")

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# calculate % unique per round
unique <- 100-diversity

# create data frame for plotting
unique_rounds <- data.frame(
  roundNums = rounds,
  uniquePercent = unique
)

# line plot with % of unique sequences per round
ggplot(data=unique_rounds, aes(roundNums, uniquePercent)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylab("% Unique Sequences") +
  xlab("Round Number")

# line plots with round 8 clustering info
# counts of sequences
ggplot(data=peaks_identity, aes(x = cluster_ID, y = sequence_counts)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_x_continuous(breaks = breaks_pretty()) +
  ylab("Raw Counts") +
  xlab("Cluster Number")

# counts of reads
ggplot(data=peaks_identity, aes(x = cluster_ID, y = read_counts)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_x_continuous(breaks = breaks_pretty()) +
  ylab("Raw Counts") +
  xlab("Cluster Number")

# avg identity
ggplot(data=peaks_identity, aes(x = cluster_ID, y = avg_ID)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_x_continuous(breaks = breaks_pretty()) +
  ylab("Average % Identity") +
  xlab("Cluster Number")

# identity info
ggplot(data=peaks_identity, aes(x = cluster_ID)) +
  geom_line(aes(y = avg_ID, color = "Avg")) +
  #geom_line(aes(y = highest_ID, color = "Highest")) +
  geom_line(aes(y = lowest_ID, color = "Lowest")) +
  geom_point(aes(y = avg_ID, color = "Avg")) +
  #geom_point(aes(y = highest_ID, color = "Highest")) +
  geom_point(aes(y = lowest_ID, color = "Lowest")) +
  scale_color_manual(name = "Counts", values = c("Avg" = safe_colors[5], "Highest" = safe_colors[9], "Lowest" = safe_colors[4])) +
  theme_bw() +
  scale_x_continuous(breaks = breaks_pretty()) +
  ylab("% Identity") +
  xlab("Cluster Number")

# line plot with fraction abundance per round for each of the top 10 families


# heatmaps with the counts (e.g., log counts) for each of the top 10 sequences per round




