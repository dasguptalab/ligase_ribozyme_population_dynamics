#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
#library(plyr)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/08_top10_sequence_abundances"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")[1:8]

# round numbers
rounds <- c(1, 2, 3, 4, 5, 6, 7, 8)

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in sequence count data
seqs_counts <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11_quantified_top10/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))

# reverse complement the sequences
#seqs_counts$sequence <- rev(chartr("ATGC","TACG",seqs_counts$sequence))

# setup data length
data_length <- 8*8
  
# data frame for the summed counts
sum_counts <- data.frame(
  run_name = rep(NA, data_length),
  counts_run_name = rep(NA, data_length),
  sum_counts = rep(NA, data_length),
  frac_abundance = rep(NA, data_length),
  run_color = rep(NA, data_length)
)

# initialize index counter
index <- 1

# loop over each run to sum top 10 sequence counts for each run and add the fraction abundances
for (run_num in 1:8) {
  # loop over each run
  for (counts_run_num in 1:8) {
    # add run name
    sum_counts$run_name[index] <- run_num
    # add counts run name
    sum_counts$counts_run_name[index] <- counts_run_num 
    # add counts
    sum_counts$sum_counts[index] <- sum(seqs_counts[seqs_counts$run_name == run_num & seqs_counts$counts_run_name == counts_run_num, "counts"])
    # add frac abundance
    sum_counts$frac_abundance[index] <- sum_counts$sum_counts[index]/quality[counts_run_num]
    # add run plotting color
    sum_counts$run_color[index] <- safe_colors[run_num]
    # update index counter
    index <- index + 1
  }
}

# add log values
sum_counts$log_sum_counts <- log(sum_counts$sum_counts)
sum_counts$log_frac_abundance <- log(sum_counts$frac_abundance)

# line plot with counts of the top 10 sequences per round for each of the rounds
cluster_abundances_plot <- ggplot(data=sum_counts, aes(x=as.character(counts_run_name), y=sum_counts, group=run_name, color=run_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Round", labels = rounds, breaks = safe_colors, guide = "legend") +
  ylab("Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/top10_sequence_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# line plot with counts of the top 10 sequences per round for each of the rounds
cluster_abundances_plot <- ggplot(data=sum_counts, aes(x=as.character(counts_run_name), y=log(sum_counts), group=run_name, color=run_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Round", labels = rounds, breaks = safe_colors, guide = "legend") +
  ylab("Log Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/top10_sequence_log_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# line plot with counts of the top 10 sequences per round for each of the rounds
cluster_abundances_plot <- ggplot(data=sum_counts, aes(x=as.character(counts_run_name), y=log(sum_counts), group=run_name, color=run_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Round", labels = rounds, breaks = safe_colors, guide = "legend") +
  ylab("Log Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/top10_sequence_log_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# line plot with counts of the top 10 sequences per round for each of the rounds
cluster_abundances_plot <- ggplot(data=sum_counts, aes(x=as.character(counts_run_name), y=frac_abundance, group=run_name, color=run_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Round", labels = rounds, breaks = safe_colors, guide = "legend") +
  ylab("Fraction Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/top10_sequence_fraction_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# line plot with counts of the top 10 sequences per round for each of the rounds
cluster_abundances_plot <- ggplot(data=sum_counts, aes(x=as.character(counts_run_name), y=100*frac_abundance, group=run_name, color=run_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Round", labels = rounds, breaks = safe_colors, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/top10_sequence_percent_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# export plotting data
write.csv(sum_counts, file = paste(out_dir, "/top10_sequence_abundances.csv", sep = ""), row.names = FALSE, quote = FALSE)
