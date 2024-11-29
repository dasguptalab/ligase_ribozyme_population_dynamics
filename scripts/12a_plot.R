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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/07a_clustered"

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# round numbers
rounds <- c(1, 2, 3, 4, 5, 6, 7, 8)

# % diversity per round
#diversity_doped <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20, 97.30, 92.40, 86.43)
diversity <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20)

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in cluster identity data
r8_peaks_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/r8_S8_L001_formatted_above9_cluster_peaks_identity_table.csv")
peaks_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/peaks_identity_table.csv")
#seqs_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/sequences_identity_table.csv")

# TO-DO: double check for duplicate data... (low ID cluster peak entries)
# read in cluster family sequence data
r8_seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tmp/10_cluster_families/07a_clustered/r8_S8_L001_cluster_families_table.csv")
#r8_peaks_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tmp/10_cluster_families/07a_clustered/r8_cluster_peaks_families_table.csv")
r8_peaks_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tmp/10_cluster_families/07a_clustered/r8_S8_L001_cluster_peaks_families_table.csv")

# read in sequence count data
seqs_counts <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_quantified_top10/07a_clustered/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))

# reverse complement the sequences
seqs_counts$sequence <- rev(chartr("ATGC","TACG",seqs_counts$sequence))

# add the fraction abundances
seqs_counts$frac_abundance <- NA
for (run_num in 1:8) {
  seqs_counts[seqs_counts$counts_run_name == run_num, "frac_abundance"] <- seqs_counts[seqs_counts$counts_run_name == run_num, "counts"]/quality[run_num]
}

# add log values
seqs_counts$log_counts <- log(seqs_counts$counts)
seqs_counts$log_frac_abundance <- log(seqs_counts$frac_abundance)

# set infinite and NA values equal to zero
is.na(seqs_counts)<-sapply(seqs_counts, is.infinite)

# re-order the data for plotting
seqs_counts <- seqs_counts[order(seqs_counts$log_counts, decreasing=TRUE),]

# setup midpoint values for plotting
seqs_counts_noNA <- seqs_counts
seqs_counts_noNA[is.na(seqs_counts_noNA)] <- 0
mid_log_counts <- max(seqs_counts_noNA$log_counts)/2
mid_log_frac_abundance <- log(min(seqs_counts_noNA[seqs_counts_noNA$frac_abundance != 0,"frac_abundance"]))/2

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
  #frac_abundance_log = rep(NA, data_length)
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
    # add log fraction abundance
    #cluster_abundances$frac_abundance_log[index] <- log(sum(r8_seqs_family[r8_seqs_family$cluster_ID == cluster_num, curr_col]))-log(quality[1])
  }
}

# add log frac abundances
cluster_abundances$log_frac_abundance <- log(cluster_abundances$frac_abundance)

# calculate % unique per round
unique_reads <- diversity

# create data frame for plotting
unique_rounds <- data.frame(
  round_nums = rounds,
  #round_nums = as.character(rounds),
  unique_percent = unique_reads
)

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

# line plot with % of unique sequences per round
unique_rounds_plot <- ggplot(data=unique_rounds, aes(round_nums, unique_percent)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylab("Percent Unique Sequences") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/percent_unique.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(unique_rounds_plot)
dev.off()

# line plots with round 8 clustering identity info
r8_peaks_identity_plot <- ggplot(data=r8_peaks_identity, aes(x = cluster_ID)) +
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
# save the plot
exportFile <- paste(out_dir, "/r8_cluster_identity.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(r8_peaks_identity_plot)
dev.off()

# line plot with fraction abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=frac_abundance, group=cluster_ID, color=cluster_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Family", labels = cluster_list_out, breaks = safe_colors, guide = "legend") +
  ylab("Fraction Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/r8_cluster_fraction_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# line plot with log of the decimal from the fraction abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=log_frac_abundance, group=cluster_ID, color=cluster_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Family", labels = cluster_list_out, breaks = safe_colors, guide = "legend") +
  ylab("Log Fraction Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/r8_cluster_log_fraction_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# export plotting data
write.csv(unique_rounds, file = paste(out_dir, "/data/percent_unique.csv", sep = ""))
write.csv(r8_peaks_identity, file = paste(out_dir, "/data/r8_cluster_identity.csv", sep = ""))
write.csv(cluster_abundances, file = paste(out_dir, "/data/r8_cluster_fraction_abundances.csv", sep = ""))
