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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/07_family_overhang_conservation"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# read in cluster sequence family data
seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/r8_S8_L001_formatted_above9_cluster_sequences_identity_table_atLeast90.csv")

# list of cluster IDs
cluster_list <- unique(seqs_family$cluster_ID)

# initialize data frame for base counts
base_counts <- data.frame(
  cluster_ID = rep(NA, 40*4),
  fam_ID = rep(NA, 40*4),
  base_ID = rep(NA, 40*4),
  base = rep(NA, 40*4),
  conservation = rep(NA, 40*4),
  comparison = rep(NA, 40*4)
)
base_counts_out <- data.frame()

# loop over each cluster
for (cluster_num in min(cluster_list):max(cluster_list)) {
  # convert list of sequences into a matrix
  seqs_matrix <- do.call(rbind, type.convert(strsplit(seqs_family[seqs_family$cluster_ID == cluster_num, "sequence"], ""), as.is = TRUE))
  # loop over each base
  for (base_num in 1:40) {
    # update indicies
    index <- ((base_num-1)*4)+1
    index_max <- index+3
    # set family number
    fam_num <- cluster_num+1
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
  # subset the sequence family data to the overhang bases
  base_counts_subset <- base_counts[base_counts$base_ID >= 16 & base_counts$base_ID <= 25,]
  # set round plot title
  run_title <- paste("Family", fam_num, "Overhang Base Conservation")
  # create heatmap of base conservation
  base_counts_plot <- ggplot(data = base_counts_subset, aes(reorder(as.character(base_ID), base_ID), base, fill= conservation_na)) + 
    theme_bw() +
    geom_tile(colour = "black") +
    #expected_overhang <- c("U","G","C","G","G","A","A","U","G","C")
    annotate("rect", 
             xmin = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5), xmax = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), 
             ymin = c(3.5, 2.5, 1.5, 2.5, 2.5, 0.5, 0.5, 3.5, 2.5, 1.5), ymax = c(4.5, 3.5, 2.5, 3.5, 3.5, 1.5, 1.5, 4.5, 3.5, 2.5), 
             colour = safe_colors[6], 
             fill = "transparent", 
             linewidth = 1) +
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
                         na.value = "white") +
    geom_text(aes(label = round(conservation_na, digits = 2)), color = "white", size = 3)
  # save the plot
  exportFile <- paste(out_dir, "/family_overhang_percent_abundance", fam_num, ".png", sep = "")
  png(exportFile, units="in", width=5, height=5, res=300)
  print(base_counts_plot)
  dev.off()
  # add current family data to outputs
  base_counts_out <- rbind(base_counts_out, base_counts_subset)
}

# export plotting data
write.csv(base_counts_out, file = paste(out_dir, "/family_overhang_conservation.csv", sep = ""), row.names = FALSE, quote = FALSE)
