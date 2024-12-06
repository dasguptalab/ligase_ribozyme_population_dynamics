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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/03_family_counts"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

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

# read in cluster peaks data
r8_peaks <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized/07a_clustered/r8_S8_L001_formatted_above9_cluster_peaks_table.csv")

# read in cluster identity data
r8_seqs_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/r8_S8_L001_formatted_above9_cluster_sequences_identity_table_atLeast90.csv")

# read in cluster family sequence data
r8_seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_families/r8_S8_L001_counts_plot_table.csv")

# retain peak counts for sequences with ID >= 90%
r8_peaks_family <- merge(r8_peaks, r8_seqs_family)

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
    #curr_col <- paste("r", run_num, "_counts", sep="")
    # add fraction abundance
    #cluster_abundances$frac_abundance[index] <- sum(r8_seqs_family[r8_seqs_family$cluster_ID == cluster_num, curr_col])/quality[run_num]
    cluster_abundances$frac_abundance[index] <- sum(r8_seqs_family[r8_seqs_family$cluster_ID == cluster_num & r8_seqs_family$counts_run_name == run_num, "counts"])/quality[run_num]
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
  fam_num = rep(NA, counts_length),
  counts = rep(NA, counts_length)
)

# data frame for peak counts
peak_counts <- data.frame()

# seq numbers for plotting
for (round_index in 1:8) {
  # set the index
  index <- ((round_index-1)*8)+1
  # setup the column name
  #curr_col <- paste("r", round_index, "_counts", sep="")
  # add round num
  subset_peak_counts$round_num <- rep(round_index, 13)
  # add seq num
  subset_peak_counts$fam_num <- seq(from = 1, to = 13, by = 1)
  # add peak seq counts
  subset_peak_counts$counts <- r8_peaks_family[r8_peaks_family$counts_run_name == round_index, "counts"]
  # add rows
  peak_counts <- rbind(peak_counts, subset_peak_counts)
}

# add the fraction abundances
peak_counts$frac_abundance <- NA
for (run_num in 1:8) {
  peak_counts[peak_counts$round_num == run_num, "frac_abundance"] <- peak_counts[peak_counts$round_num == run_num, "counts"]/quality[run_num]
}

# change zeros to NAs for plotting
peak_counts$frac_abundance_na <- ifelse(peak_counts$frac_abundance == 0, NA, peak_counts$frac_abundance)

# heatmaps with the log counts for each of sequence families from round 8
peak_counts_plot <- ggplot(data = peak_counts, aes(as.character(round_num), reorder(as.character(fam_num), fam_num, decreasing = TRUE), fill= log(counts))) + 
  theme_bw() +
  geom_tile(colour = "black") +
  ylab("Family ID") +
  xlab("Round Number") +
  scale_fill_gradient2(name = "Log Counts",
                       low = safe_colors[3],
                       mid = safe_colors[4],
                       high = safe_colors[5],
                       midpoint = log(max(peak_counts$counts))/2,
                       na.value = "white")
# save the plot
exportFile <- paste(out_dir, "/family_log_counts.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(peak_counts_plot)
dev.off()

# heatmap with the fraction abundances for each of sequence families from round 8
peak_counts_plot <- ggplot(data = peak_counts, aes(as.character(round_num), reorder(as.character(fam_num), fam_num, decreasing = TRUE), fill= frac_abundance_na)) + 
  theme_bw() +
  geom_tile(colour = "black") +
  ylab("Family ID") +
  xlab("Round Number") +
  scale_fill_gradient2(name = "FA",
                       low = safe_colors[3],
                       mid = safe_colors[4],
                       high = safe_colors[5],
                       midpoint = max(peak_counts$frac_abundance)/2,
                       na.value = "white")
# save the plot
exportFile <- paste(out_dir, "/family_fraction_abundance.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(peak_counts_plot)
dev.off()

# heatmap with the fraction abundances for each of sequence families from round 8
peak_counts_plot <- ggplot(data = peak_counts, aes(as.character(round_num), reorder(as.character(fam_num), fam_num, decreasing = TRUE), fill= 100*frac_abundance_na)) + 
  theme_bw() +
  geom_tile(colour = "black") +
  ylab("Family ID") +
  xlab("Round Number") +
  scale_fill_gradient2(name = "PA",
                       low = safe_colors[3],
                       mid = safe_colors[4],
                       high = safe_colors[5],
                       midpoint = 100*(max(peak_counts$frac_abundance)/2),
                       na.value = "white")
# save the plot
exportFile <- paste(out_dir, "/family_percent_abundance.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(peak_counts_plot)
dev.off()

# export plotting data
write.csv(peak_counts, file = paste(out_dir, "/family_counts.csv", sep = ""), row.names = FALSE, quote = FALSE)
