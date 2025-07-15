#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
library(stringr)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/F4a_family_base_conservation"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# To-do: retrieve seqs data from cluster summaries?
# read in cluster sequence family data
#seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/r8_S8_L001_formatted_above9_cluster_sequences_identity_table_atLeast90.csv")

# read in sequences that have at least 90% identity to any peak
seqs_identities <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/10_family_identification/family_identities_atLeast90.csv")

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
cluster_list <- unique(seqs_family$peak_cluster_IDs)

# list of cluster IDs and abundances
cluster_read_counts <- c(40999, 256601, 34117, 82079, 14439, 26345, 6555, 12539, 15921, 2572, 4227, 8543, 4043)
cluster_seq_counts <- c(108, 1114, 98, 116, 269, 85, 33, 49, 75, 153, 275, 286, 271)
cluster_read_abun <- 100*cluster_read_counts/889374

# create cluster data frame
cluster_data <- data.frame(
  cluster_ID = cluster_list,
  read_counts = cluster_read_counts,
  seq_counts = cluster_seq_counts,
  read_abun = cluster_read_abun
)

# sort cluster data
cluster_data <- cluster_data[order(cluster_data$read_abun, decreasing = TRUE),]  

# add family numbers
cluster_data$fam_num <- seq(from = 1, to = 13, by = 1)

# list of cluster IDs
cluster_list <- unique(seqs_family$cluster_ID)

# initialize data frame for base counts
base_counts <- data.frame(
  cluster_ID = rep(NA, 40*4),
  fam_ID = rep(NA, 40*4),
  base_ID = rep(NA, 40*4),
  base = rep(NA, 40*4),
  conservation = rep(NA, 40*4)
)
base_counts_out <- data.frame()

# initialize data frame for sequences
seqs_matrix <- data.frame()

# loop over each cluster
for (cluster_num in min(cluster_list):max(cluster_list)) {
#for (cluster_num in 1:1) {
  #cluster_num <- 1
  # convert list of sequences into a matrix
  seqs_matrix <- do.call(rbind, type.convert(strsplit(seqs_family[seqs_family$cluster_ID == cluster_num, "sequence"], ""), as.is = TRUE))
  # set family number
  fam_num <- cluster_data[cluster_data$cluster_ID == cluster_num, "fam_ID"]
  # loop over each base
  for (base_num in 1:40) {
    # update indicies
    index <- ((base_num-1)*4)+1
    index_max <- index+3
    # add cluster ID
    base_counts$cluster_ID[index:index_max] <- rep(cluster_num, 4)
    # add family ID
    base_counts$fam_ID[index:index_max] <- rep(fam_num, 4)
    # add base ID
    base_counts$base_ID[index:index_max] <- rep(base_num, 4)
    # add percent conservation of each base character
    base_counts$conservation[index] <- 100*sum(str_count(seqs_matrix[,base_num], "A"))/nrow(seqs_matrix)
    base_counts$conservation[index+1] <- 100*sum(str_count(seqs_matrix[,base_num], "C"))/nrow(seqs_matrix)
    base_counts$conservation[index+2] <- 100*sum(str_count(seqs_matrix[,base_num], "G"))/nrow(seqs_matrix)
    base_counts$conservation[index_max] <- 100*sum(str_count(seqs_matrix[,base_num], "T"))/nrow(seqs_matrix)
    # add each base character
    base_counts$base[index] <- "A"
    base_counts$base[index+1] <- "C"
    base_counts$base[index+2] <- "G"
    base_counts$base[index_max] <- "U"
  }
  # change zeros to NAs for plotting
  base_counts$conservation_na <- ifelse(base_counts$conservation == 0, NA, base_counts$conservation)
  # change the bases to factors for plotting
  #base_counts$base <- factor(base_counts$base, levels=c("U", "G", "C", "A"))
  # set family plot title
  run_title <- paste("Family", fam_num, "Base Conservation")
  # create heatmap of base conservation
  base_counts_plot <- ggplot(data = base_counts, aes(reorder(as.character(base_ID), base_ID), base, fill= conservation_na)) + 
    theme_classic() +
    geom_tile(colour = "black") +
    # left P2
    annotate("rect", xmin = c(5.5), xmax = c(12.5), ymin = c(0.5), ymax = c(4.5), 
             colour = safe_colors[2], fill = "transparent", linewidth = 1.25) +
    # right P2
    annotate("rect", xmin = c(25.5), xmax = c(32.5), ymin = c(0.5), ymax = c(4.5), 
             colour = safe_colors[2], fill = "transparent", linewidth = 1.25) +
    # overhang compliment
    annotate("rect", xmin = c(15.5), xmax = c(23.5), ymin = c(0.5), ymax = c(4.5),
             #xmin = c(15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5), xmax = c(16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5), 
             #ymin = c(3.5, 2.5, 1.5, 2.5, 2.5, 0.5, 0.5, 3.5, 2.5, 1.5), ymax = c(4.5, 3.5, 2.5, 3.5, 3.5, 1.5, 1.5, 4.5, 3.5, 2.5), 
             #xmin = c(15.5, 16.5, 17.5, 18.5, 19.5), xmax = c(16.5, 17.5, 18.5, 19.5, 20.5), 
             #ymin = c(3.5, 2.5, 1.5, 2.5, 2.5), ymax = c(4.5, 3.5, 2.5, 3.5, 3.5), 
             colour = safe_colors[1], fill = "transparent", linewidth = 1.25) +
    ylab("Base") +
    xlab("Base Number") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), 
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 18)) +
    ggtitle(run_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 18)) +
    scale_fill_gradient2(name = "PA",
                         low = safe_colors[3],
                         mid = safe_colors[4],
                         high = safe_colors[5],
                         midpoint = max(base_counts$conservation)/2,
                         na.value = "white") +
    coord_fixed() +
    geom_text(aes(label = round(conservation_na, digits = 2)), color = "white", size = 4)
  # save the plot
  exportFile <- paste(out_dir, "/family", fam_num, "_base_conservation.png", sep = "")
  png(exportFile, units="in", width=20, height=3, res=300)
  print(base_counts_plot)
  dev.off()
  # add current family data to outputs
  base_counts_out <- rbind(base_counts_out, base_counts)
}

# export plotting data
write.csv(base_counts_out, file = paste(out_dir, "/family_base_conservation.csv", sep = ""), row.names = FALSE, quote = FALSE)
