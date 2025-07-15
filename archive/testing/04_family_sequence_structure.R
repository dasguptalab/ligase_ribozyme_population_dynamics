#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# install logo tool
# https://omarwagih.github.io/ggseqlogo/
#devtools::install_github("omarwagih/ggseqlogo")

# import libraries
library(ggplot2)
#library(scales)
library(rcartocolor)
library(stringr)
require(ggseqlogo)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/04_family_sequence_structure"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# read in cluster sequence family data
seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/r8_S8_L001_formatted_above9_cluster_sequences_identity_table_atLeast90.csv")

# convert Ts to Us
seqs_input$sequence <- gsub("T", "U", seqs_input$sequence)

# read in sequences that have at least 90% identity to any peak
seqs_identities <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/00a_family_identification_above9/family_identities_above9_atLeast90.csv")

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
fam_list_out <- seq(1, 13)

# list of cluster IDs in order of abundance in round 8
r8_fams <- data.frame(
  fam_ID = fam_list_out,
  cluster_ID = c(1, 3, 0, 2, 5, 8, 4, 7, 11, 6, 10, 12, 9)
)

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

# set number of bases
num_bases <- 6

# set number of base pairs
num_base_pairs <- 7

# set number of base pairings
num_pairs <- 16

# set data length
data_length <- num_base_pairs*num_pairs

# initialize data frames and variables for base pairing counts
base_pair_counts <- data.frame(
  cluster_ID = rep(NA, data_length),
  fam_ID = rep(NA, data_length),
  base_ID = rep(NA, data_length),
  base = rep(NA, data_length),
  counts = rep(NA, data_length),
  conservation = rep(NA, data_length)
)
base_pair_counts_out <- data.frame()

# loop over each cluster
for (cluster_num in min(cluster_list):max(cluster_list)) {
  # convert list of sequences into a matrix
  seqs_matrix <- do.call(rbind, type.convert(strsplit(seqs_family[seqs_family$cluster_ID == cluster_num, "sequence"], ""), as.is = TRUE))
  # set family number
  fam_num <- r8_fams[r8_fams$cluster_ID == cluster_num, "fam_ID"]
  # initialize right stem base counter
  base_num_left <- 6
  base_num_right <- 32
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
    base_counts$conservation[index_max] <- 100*sum(str_count(seqs_matrix[,base_num], "U"))/nrow(seqs_matrix)
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
    # add current family data to outputs
    base_counts_out <- rbind(base_counts_out, base_counts)
  # loop over each stem or overhang base
  for (base_num in 1:num_base_pairs) {
    # update indexing variables
    index <- ((base_num-1)*num_pairs)+1
    index_max <- index+(num_pairs-1)
    # add cluster ID
    base_pair_counts$cluster_ID[index:index_max] <- rep(cluster_num, num_pairs)
    # add family ID
    base_pair_counts$fam_ID[index:index_max] <- rep(fam_num, num_pairs)
    # add base ID
    base_pair_counts$base_ID[index:index_max] <- rep(paste(base_num_left, base_num_right, sep="_"), num_pairs)
    # add counts of each base pairing
    base_pair_counts$counts[index] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "A" & seqs_matrix[,base_num_right] == "U", base_num_left]))
    base_pair_counts$counts[index+1] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "U" & seqs_matrix[,base_num_right] == "A", base_num_left]))
    base_pair_counts$counts[index+2] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "G" & seqs_matrix[,base_num_right] == "C", base_num_left]))
    base_pair_counts$counts[index+3] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "C" & seqs_matrix[,base_num_right] == "G", base_num_left]))
    base_pair_counts$counts[index+4] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "U" & seqs_matrix[,base_num_right] == "G", base_num_left]))
    base_pair_counts$counts[index+5] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "G" & seqs_matrix[,base_num_right] == "U", base_num_left]))
    base_pair_counts$counts[index+6] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "U" & seqs_matrix[,base_num_right] == "C", base_num_left]))
    base_pair_counts$counts[index+7] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "C" & seqs_matrix[,base_num_right] == "U", base_num_left]))
    base_pair_counts$counts[index+8] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "A" & seqs_matrix[,base_num_right] == "C", base_num_left]))
    base_pair_counts$counts[index+9] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "C" & seqs_matrix[,base_num_right] == "A", base_num_left]))
    base_pair_counts$counts[index+10] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "A" & seqs_matrix[,base_num_right] == "G", base_num_left]))
    base_pair_counts$counts[index+11] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "G" & seqs_matrix[,base_num_right] == "A", base_num_left]))
    base_pair_counts$counts[index+12] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "A" & seqs_matrix[,base_num_right] == "A", base_num_left]))
    base_pair_counts$counts[index+13] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "C" & seqs_matrix[,base_num_right] == "C", base_num_left]))
    base_pair_counts$counts[index+14] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "G" & seqs_matrix[,base_num_right] == "G", base_num_left]))
    base_pair_counts$counts[index_max] <- sum(str_count(seqs_matrix[seqs_matrix[,base_num_left] == "U" & seqs_matrix[,base_num_right] == "U", base_num_left]))
    # add conservation of each base pairing
    base_pair_counts$conservation[index] <- 100*base_pair_counts$counts[index]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+1] <- 100*base_pair_counts$counts[index+1]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+2] <- 100*base_pair_counts$counts[index+2]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+3] <- 100*base_pair_counts$counts[index+3]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+4] <- 100*base_pair_counts$counts[index+4]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+5] <- 100*base_pair_counts$counts[index+5]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+6] <- 100*base_pair_counts$counts[index+6]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+7] <- 100*base_pair_counts$counts[index+7]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+8] <- 100*base_pair_counts$counts[index+8]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+9] <- 100*base_pair_counts$counts[index+9]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+10] <- 100*base_pair_counts$counts[index+10]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+11] <- 100*base_pair_counts$counts[index+11]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+12] <- 100*base_pair_counts$counts[index+12]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+13] <- 100*base_pair_counts$counts[index+13]/nrow(seqs_matrix)
    base_pair_counts$conservation[index+14] <- 100*base_pair_counts$counts[index+14]/nrow(seqs_matrix)
    base_pair_counts$conservation[index_max] <- 100*base_pair_counts$counts[index_max]/nrow(seqs_matrix)
    # add each base pairing
    base_pair_counts$base[index] <- "A-U"
    base_pair_counts$base[index+1] <- "U-A"
    base_pair_counts$base[index+2] <- "G-C"
    base_pair_counts$base[index+3] <- "C-G"
    base_pair_counts$base[index+4] <- "U-G"
    base_pair_counts$base[index+5] <- "G-U"
    base_pair_counts$base[index+6] <- "U-C"
    base_pair_counts$base[index+7] <- "C-U"
    base_pair_counts$base[index+8] <- "A-C"
    base_pair_counts$base[index+9] <- "C-A"
    base_pair_counts$base[index+10] <- "A-G"
    base_pair_counts$base[index+11] <- "G-A"
    base_pair_counts$base[index+12] <- "A-A"
    base_pair_counts$base[index+13] <- "C-C"
    base_pair_counts$base[index+14] <- "G-G"
    base_pair_counts$base[index_max] <- "U-U"
    # increment left stem base counter
    base_num_left <- base_num_left+1
    # decrement right stem base counter
    base_num_right <- base_num_right-1
  }
  # change zeros to NAs for plotting
  base_pair_counts$conservation_na <- ifelse(base_pair_counts$conservation == 0, NA, base_pair_counts$conservation)
  # change the base pair numbers to factors for plotting
  base_pair_counts$base_ID <- factor(base_pair_counts$base_ID, levels=c("6_32", "7_31", "8_30", "9_29", "10_28", "11_27", "12_26"))
  # add current family data to outputs
  base_pair_counts_out <- rbind(base_pair_counts_out, base_pair_counts)
  # set family plot title
  run_title <- paste("Family", fam_num, "Sequence Logo")
  # create sequence logo
  logo_plot <- ggplot() + 
    geom_logo(seqs_family[seqs_family$cluster_ID == cluster_num, "sequence"], method = 'prob') + 
    # left stem
    annotate('text', x=9, y=1.05, label='Left Stem') + 
    annotate('rect', xmin = 5.5, xmax = 12.5, ymin = -0.01, ymax = 1.01, 
             alpha = .1, col=safe_colors[2], fill=safe_colors[3], linewidth = 1) + 
    # right stem
    annotate('text', x=29, y=1.05, label='Right Stem') + 
    annotate('rect', xmin = 25.5, xmax = 32.5, ymin = -0.01, ymax = 1.01, 
             alpha = .1, col=safe_colors[2], fill=safe_colors[3], linewidth = 1) + 
    # overhang
    annotate('text', x=18, y=1.05, label='Overhang') + 
    annotate('rect', xmin = 15.5, xmax = 20.5, ymin = -0.01, ymax = 1.01, 
             alpha = .1, col=safe_colors[1], fill=safe_colors[3], linewidth = 1) + 
    theme_logo() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
    ggtitle(run_title) +
    theme(plot.title = element_text(hjust = 0.5))
  # set family plot title
  run_title <- paste("Family", fam_num, "Base Conservation")
  # create heatmap of base conservation
  base_counts_plot <- ggplot(data = base_counts, aes(reorder(as.character(base_ID), base_ID), base, fill= conservation_na)) + 
    theme_bw() +
    geom_tile(colour = "black") +
    # left stem
    annotate("rect", xmin = c(5.5), xmax = c(12.5), ymin = c(0.5), ymax = c(4.5), 
             colour = safe_colors[2], fill = "transparent", linewidth = 1) +
    # right stem
    annotate("rect", xmin = c(25.5), xmax = c(32.5), ymin = c(0.5), ymax = c(4.5), 
             colour = safe_colors[2], fill = "transparent", linewidth = 1) +
    # overhang
    annotate("rect", 
             xmin = c(15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5), xmax = c(16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5), 
             ymin = c(3.5, 2.5, 1.5, 2.5, 2.5, 0.5, 0.5, 3.5, 2.5, 1.5), ymax = c(4.5, 3.5, 2.5, 3.5, 3.5, 1.5, 1.5, 4.5, 3.5, 2.5), 
             colour = safe_colors[1], fill = "transparent", linewidth = 1) +
    ylab("Base") +
    xlab("Base Number") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
    ggtitle(run_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(name = "PA",
                         low = safe_colors[3],
                         mid = safe_colors[4],
                         high = safe_colors[5],
                         midpoint = max(base_counts$conservation)/2,
                         na.value = "white") #+
    #coord_fixed()
  # set round plot title
  run_title <- paste("Family", fam_num, "Stem Base Pairing Conservation")
  # create heatmap of base conservation
  base_pair_counts_plot <- ggplot(data = base_pair_counts_out, aes(reorder(base_ID, base_ID), base, fill= conservation_na)) + 
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
                         midpoint = max(base_pair_counts_out$conservation)/2,
                         na.value = "white") #+
    #coord_fixed()
    #geom_text(aes(label = round(conservation_na, digits = 2)), color = "white", size = 4)
  # combine plots
  suppressMessages( require(cowplot) )
  grid_plot <- plot_grid(logo_plot, base_counts_plot, base_pair_counts_plot,  ncol = 1, align = 'v')
  #grid_plot <- plot_grid(logo_plot, base_counts_plot,  ncol = 1, align = 'v')
  # save the plot
  exportFile <- paste(out_dir, "/family", fam_num, "_sequence_structure.png", sep = "")
  png(exportFile, units="in", width=10, height=12, res=300)
  print(grid_plot)
  dev.off()
}

# export plotting data
write.csv(seqs_family, file = paste(out_dir, "/family_sequences.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(base_counts_out, file = paste(out_dir, "/family_base_conservation.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(base_pair_counts_out, file = paste(out_dir, "/family_stem_base_pairings.csv", sep = ""), row.names = FALSE, quote = FALSE)
