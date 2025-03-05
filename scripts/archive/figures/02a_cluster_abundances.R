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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/02a_cluster_abundances"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
unique <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)
#above_9 <- c(5, 3, 5, 4, 283, 4001, 2703, 2100)

# list of cluster IDs and abundances
cluster_list <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
cluster_read_counts <- c(40999, 256601, 34117, 82079, 14439, 26345, 6555, 12539, 15921, 2572, 4227, 8543, 4043)
cluster_seq_counts <- c(108, 1114, 98, 116, 269, 85, 33, 49, 75, 153, 275, 286, 271)
cluster_read_abun <- 100*cluster_read_counts/889374

# create cluster data frame
cluster_data <- data.frame(
  cluster_ID = cluster_list,
  read_counts = cluster_read_counts,
  seq_counts = cluster_seq_counts,
  read_abun = cluster_read_abun,
  cluster_color = safe_colors
)

# sort cluster data
cluster_data <- cluster_data[order(cluster_data$read_abun, decreasing = TRUE),]  

# add family numbers
cluster_data$fam_num <- seq(from = 1, to = 13, by = 1)
  
# read in cluster family sequence data
r8_seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_families/r8_S8_L001_counts_plot_table.csv")

# setup data frame length
data_length <- (max(cluster_list)+1)*8

# data frame of cluster abundances
cluster_abundances <- data.frame(
  cluster_ID = rep(NA, data_length), 
  run_name = rep(NA, data_length),
  cluster_color = rep(NA, data_length),
  read_counts = rep(NA, data_length),
  unique_counts = rep(NA, data_length),
  fam_ID = rep(NA, data_length)
)

# calculate fraction abundance per round
for (cluster_num in 0:max(cluster_list)) {
  # loop over each run
  for (run_num in 1:8) {
    # set the index
    index <- run_num+(cluster_num*8)
    # add family number
    cluster_abundances$fam_ID[index] <- cluster_data[cluster_data$cluster_ID == cluster_num, "fam_num"]
    # add cluster ID
    cluster_abundances$cluster_ID[index] <- cluster_num
    # add run name
    cluster_abundances$run_name[index] <- run_num
    # add cluster plotting color
    cluster_abundances$cluster_color[index] <- safe_colors[cluster_num+1]
    # setup the column name
    #curr_col <- paste("r", run_num, "_counts", sep="")
    # To-do: fix unique counts... currently counts just number of clusters
    # add counts
    cluster_abundances$unique_counts[index] <- nrow(r8_seqs_family[r8_seqs_family$cluster_ID == cluster_num,])
    cluster_abundances$read_counts[index] <- sum(r8_seqs_family[r8_seqs_family$cluster_ID == cluster_num & r8_seqs_family$counts_run_name == run_num, "counts"])
    # add fraction abundances
    cluster_abundances$unique_frac_abundance[index] <- cluster_abundances$unique_counts[index]/unique[run_num]
    cluster_abundances$read_frac_abundance[index] <- cluster_abundances$read_counts[index]/quality[run_num]
  }
}

# add percentage abundances
cluster_abundances$unique_perc_abundance <- 100*cluster_abundances$unique_frac_abundance
cluster_abundances$read_perc_abundance <- 100*cluster_abundances$read_frac_abundance

# setup length of counts
counts_length <- 13*8

# sort plotting data
cluster_abundances <- cluster_abundances[order(cluster_abundances$fam_ID, decreasing = FALSE),]  

# line plot with percent abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=unique_perc_abundance, group=fam_ID, color=cluster_color))+
  geom_line(size = 1) +
  geom_point() +
  theme_classic() +
  scale_color_identity(name = "Family", labels = cluster_abundances$fam_ID, breaks = cluster_abundances$cluster_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/family_cluster_unique_percent_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(cluster_abundances_plot)
dev.off()

# line plot with percent abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=read_perc_abundance, group=fam_ID, color=cluster_color))+
  geom_line(size = 1) +
  geom_point() +
  theme_classic() +
  scale_color_identity(name = "Family", labels = cluster_abundances$fam_ID, breaks = cluster_abundances$cluster_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/family_cluster_read_percent_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(cluster_abundances_plot)
dev.off()

# export plotting data
write.csv(cluster_data, file = paste(out_dir, "/r8_cluster_count_data.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(cluster_abundances, file = paste(out_dir, "/family_cluster_abundances.csv", sep = ""), row.names = FALSE, quote = FALSE)
