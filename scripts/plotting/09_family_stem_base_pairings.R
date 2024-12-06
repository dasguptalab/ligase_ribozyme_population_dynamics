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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/09_family_stem_base_pairings"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# read in cluster sequence family data
seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/r8_S8_L001_formatted_above9_cluster_sequences_identity_table_atLeast90.csv")

# list of cluster IDs
cluster_list <- unique(seqs_family$cluster_ID)

# set number of bases
num_base_pairs <- 7

# set number of base pairings
num_pairs <- 16

# set data length
data_length <- num_base_pairs*num_pairs

# initialize data frames and variables for base pairing counts
base_counts <- data.frame(
  cluster_ID = rep(NA, data_length),
  fam_ID = rep(NA, data_length),
  base_ID = rep(NA, data_length),
  base = rep(NA, data_length),
  counts = rep(NA, data_length),
  conservation = rep(NA, data_length)
)
base_counts_out <- data.frame()

# loop over each cluster
for (cluster_num in min(cluster_list):max(cluster_list)) {
  # convert list of sequences into a matrix
  seqs_matrix <- do.call(rbind, type.convert(strsplit(seqs_family[seqs_family$cluster_ID == cluster_num, "sequence"], ""), as.is = TRUE))
  # initialize right stem base counter
  base_num_left <- 6
  base_num_right <- 32
  # loop over each stem base
  for (base_num in 1:num_base_pairs) {
    # update indexing variables
    index <- ((base_num-1)*num_pairs)+1
    index_max <- index+(num_pairs-1)
    # set family number
    fam_num <- cluster_num+1
    # add cluster ID
    base_counts$cluster_ID[index:index_max] <- rep(cluster_num, num_pairs)
    # add family ID
    base_counts$fam_ID[index:index_max] <- rep(fam_num, num_pairs)
    # add base ID
    base_counts$base_ID[index:index_max] <- rep(paste(base_num_left, base_num_right, sep="_"), num_pairs)
    # add counts of each base pairing
    base_counts$counts[index] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "A" & seqs_matrix[,base_num_right] == "T", base_num_left]))
    base_counts$counts[index+1] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "T" & seqs_matrix[,base_num_right] == "A", base_num_left]))
    base_counts$counts[index+2] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "G" & seqs_matrix[,base_num_right] == "C", base_num_left]))
    base_counts$counts[index+3] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "C" & seqs_matrix[,base_num_right] == "G", base_num_left]))
    base_counts$counts[index+4] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "T" & seqs_matrix[,base_num_right] == "G", base_num_left]))
    base_counts$counts[index+5] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "G" & seqs_matrix[,base_num_right] == "T", base_num_left]))
    base_counts$counts[index+6] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "T" & seqs_matrix[,base_num_right] == "C", base_num_left]))
    base_counts$counts[index+7] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "C" & seqs_matrix[,base_num_right] == "T", base_num_left]))
    base_counts$counts[index+8] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "A" & seqs_matrix[,base_num_right] == "C", base_num_left]))
    base_counts$counts[index+9] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "C" & seqs_matrix[,base_num_right] == "A", base_num_left]))
    base_counts$counts[index+10] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "A" & seqs_matrix[,base_num_right] == "G", base_num_left]))
    base_counts$counts[index+11] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "G" & seqs_matrix[,base_num_right] == "A", base_num_left]))
    base_counts$counts[index+12] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "A" & seqs_matrix[,base_num_right] == "A", base_num_left]))
    base_counts$counts[index+13] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "C" & seqs_matrix[,base_num_right] == "C", base_num_left]))
    base_counts$counts[index+14] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "G" & seqs_matrix[,base_num_right] == "G", base_num_left]))
    base_counts$counts[index_max] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "T" & seqs_matrix[,base_num_right] == "T", base_num_left]))
    # add conservation of each base pairing
    base_counts$conservation[index] <- 100*base_counts$counts[index]/nrow(seqs_matrix)
    base_counts$conservation[index+1] <- 100*base_counts$counts[index+1]/nrow(seqs_matrix)
    base_counts$conservation[index+2] <- 100*base_counts$counts[index+2]/nrow(seqs_matrix)
    base_counts$conservation[index+3] <- 100*base_counts$counts[index+3]/nrow(seqs_matrix)
    base_counts$conservation[index+4] <- 100*base_counts$counts[index+4]/nrow(seqs_matrix)
    base_counts$conservation[index+5] <- 100*base_counts$counts[index+5]/nrow(seqs_matrix)
    base_counts$conservation[index+6] <- 100*base_counts$counts[index+6]/nrow(seqs_matrix)
    base_counts$conservation[index+7] <- 100*base_counts$counts[index+7]/nrow(seqs_matrix)
    base_counts$conservation[index+8] <- 100*base_counts$counts[index+8]/nrow(seqs_matrix)
    base_counts$conservation[index+9] <- 100*base_counts$counts[index+9]/nrow(seqs_matrix)
    base_counts$conservation[index+10] <- 100*base_counts$counts[index+10]/nrow(seqs_matrix)
    base_counts$conservation[index+11] <- 100*base_counts$counts[index+11]/nrow(seqs_matrix)
    base_counts$conservation[index+12] <- 100*base_counts$counts[index+12]/nrow(seqs_matrix)
    base_counts$conservation[index+13] <- 100*base_counts$counts[index+13]/nrow(seqs_matrix)
    base_counts$conservation[index+14] <- 100*base_counts$counts[index+14]/nrow(seqs_matrix)
    base_counts$conservation[index_max] <- 100*base_counts$counts[index_max]/nrow(seqs_matrix)
    # add each base pairing
    base_counts$base[index] <- "A-U"
    base_counts$base[index+1] <- "U-A"
    base_counts$base[index+2] <- "G-C"
    base_counts$base[index+3] <- "C-G"
    base_counts$base[index+4] <- "U-G"
    base_counts$base[index+5] <- "G-U"
    base_counts$base[index+6] <- "U-C"
    base_counts$base[index+7] <- "C-U"
    base_counts$base[index+8] <- "A-C"
    base_counts$base[index+9] <- "C-A"
    base_counts$base[index+10] <- "A-G"
    base_counts$base[index+11] <- "G-A"
    base_counts$base[index+12] <- "A-A"
    base_counts$base[index+13] <- "C-C"
    base_counts$base[index+14] <- "G-G"
    base_counts$base[index_max] <- "U-U"
    # increment left stem base counter
    base_num_left <- base_num_left+1
    # decrement right stem base counter
    base_num_right <- base_num_right-1
  }
  # change zeros to NAs for plotting
  base_counts$conservation_na <- ifelse(base_counts$conservation == 0, NA, base_counts$conservation)
  # change the base pair numbers to factors for plotting
  base_counts$base_ID <- factor(base_counts$base_ID, levels=c("6_32", "7_31", "8_30", "9_29", "10_28", "11_27", "12_26"))
  # set round plot title
  run_title <- paste("Family", fam_num, "Stem Base Pairing Conservation")
  # create heatmap of base conservation
  base_counts_plot <- ggplot(data = base_counts, aes(reorder(base_ID, base_ID), base, fill= conservation_na)) + 
    theme_bw() +
    geom_tile(colour = "black") +
    ylab("Base Pairings") +
    xlab("Stem Base Numbers") +
    annotate("rect", 
             # base pairing rules: A with T (U) & G with C
             xmin = c(0.5, 0.5, 0.5, 0.5), xmax = c(7.5, 7.5, 7.5, 7.5), 
             ymin = c(3.5, 12.5, 6.5, 9.5), ymax = c(4.5, 13.5, 7.5, 10.5), 
             colour = safe_colors[2], fill = "transparent", linewidth = 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
    ggtitle(run_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(name = "PA",
                         low = safe_colors[3],
                         mid = safe_colors[4],
                         high = safe_colors[5],
                         midpoint = max(base_counts$conservation)/2,
                         na.value = "white") +
  geom_text(aes(label = round(conservation_na, digits = 2)), color = "white", size = 4)
  # save the plot
  exportFile <- paste(out_dir, "/family", fam_num, "_stem_base_pairings.png", sep = "")
  png(exportFile, units="in", width=5, height=5, res=300)
  print(base_counts_plot)
  dev.off()
  # add current family data to outputs
  base_counts_out <- rbind(base_counts_out, base_counts)
}

# export plotting data
write.csv(base_counts_out, file = paste(out_dir, "/family_stem_base_pairings.csv", sep = ""), row.names = FALSE, quote = FALSE)
