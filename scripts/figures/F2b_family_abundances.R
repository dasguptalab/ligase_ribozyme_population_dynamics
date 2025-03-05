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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/F2b_family_abundances"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in cluster family sequence data
#seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_families/r8_S8_L001_counts_plot_table.csv")
#seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/r8_S8_L001_formatted_above9_cluster_sequences_identity_table_atLeast90.csv")
seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))

# read in sequences that have at least 90% identity to any peak
seqs_identities <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/10a_family_identification/family_identities_atLeast90.csv")

# subset to keep round 8 data
#seqs_identities <- seqs_identities[seqs_identities$run_name == 8,]

# initialize data frame
seqs_family <- data.frame()
peak_cluster_IDs <- NULL
counts_run_ID <- NULL

# loop over each sequences that has at least 90% identity to any peak
for (seq_num in 1:nrow(seqs_identities)) {
  # keep count data for sequences that have at least 90% identity to any peak
  seqs_90_data <- seqs_input[seqs_input$sequence_ID == seqs_identities$sequence_ID[seq_num],]
  # add counts
  seqs_family <- rbind(seqs_family, seqs_90_data)
  # add peak cluster ID to vector
  peak_cluster_IDs <- c(peak_cluster_IDs, rep(seqs_identities$peak_cluster_ID[seq_num], nrow(seqs_90_data)))
  # add counts run name
  counts_run_ID <- c(counts_run_ID, rep(seqs_identities$run_name[seq_num], nrow(seqs_90_data)))
}

# add peak cluster IDs
seqs_family <- cbind(seqs_family, counts_run_ID, peak_cluster_IDs)

# list of cluster IDs
#cluster_list <- unique(seqs_family$peak_cluster_IDs)

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
#cluster_data <- cluster_data[order(cluster_data$read_abun, decreasing = TRUE),]  

# add family numbers
#cluster_data$fam_num <- seq(from = 1, to = 13, by = 1)

# setup data frame length
data_length <- (max(cluster_list)+1)*8

# data frame of cluster abundances
cluster_abundances <- data.frame(
  #fam_ID = rep(NA, data_length),
  cluster_ID = rep(NA, data_length), 
  run_name = rep(NA, data_length),
  cluster_color = rep(NA, data_length),
  counts = rep(NA, data_length),
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
    #cluster_out <- cluster_data[cluster_data$cluster_ID == cluster_num, "fam_num"]
    # add family number
    #cluster_abundances$fam_ID[index] <- paste("family", cluster_out, sep="_")
    #cluster_abundances$fam_ID[index] <- cluster_out
    # add cluster ID
    cluster_abundances$cluster_ID[index] <- cluster_num
    # add run name
    cluster_abundances$run_name[index] <- run_num
    # add cluster plotting color
    cluster_abundances$cluster_color[index] <- cluster_data[cluster_data$cluster_ID == cluster_num, "cluster_color"]
    # setup the column name
    #curr_col <- paste("r", run_num, "_counts", sep="")
    # add counts
    cluster_abundances$counts[index] <- sum(seqs_family[seqs_family$peak_cluster_ID == cluster_num & seqs_family$counts_run_ID == run_num, "counts"])
    # add fraction abundance
    cluster_abundances$frac_abundance[index] <- cluster_abundances$counts[index]/quality[run_num]
  }
}

# add percentage abundances
cluster_abundances$perc_abundance <- 100*cluster_abundances$frac_abundance

# setup length of counts
counts_length <- 13*8

# subset cluster data
cluster_abundances_subset <- cluster_abundances[cluster_abundances$run_name == 8,]
# sort cluster data
cluster_abundances_subset <- cluster_abundances_subset[order(cluster_abundances_subset$perc_abundance, decreasing = TRUE),]  
# add family numbers
cluster_abundances_subset$fam_num <- seq(from = 1, to = 13, by = 1)

# loop over each run
for (cluster_num in 0:12) {
  # update family numbers for each run
  cluster_abundances[cluster_abundances$cluster_ID == cluster_num, "fam_ID"] <- cluster_abundances_subset[cluster_abundances_subset$cluster_ID == cluster_num, "fam_num"]
}

# sort plotting data
cluster_abundances <- cluster_abundances[order(cluster_abundances$fam_ID, decreasing = FALSE),]  

# line plot with percent abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=perc_abundance, group=fam_ID, color=cluster_color))+
  geom_line(size = 1) +
  geom_point() +
  theme_classic() +
  scale_color_identity(name = "Family", labels = cluster_abundances$fam_ID, breaks = cluster_abundances$cluster_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/family_percent_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(cluster_abundances_plot)
dev.off()

# export plotting data
write.csv(cluster_data, file = paste(out_dir, "/r8_cluster_count_data.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(cluster_abundances, file = paste(out_dir, "/family_abundances.csv", sep = ""), row.names = FALSE, quote = FALSE)
