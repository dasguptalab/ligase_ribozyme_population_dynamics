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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/07_above9_overhang_conservation"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# set expected overhang sequence
#expected_overhang <- c("T","G","C","G","G","A","A","T","G","C")
expected_overhang <- c("TGCGGAATGC")

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in sequence data
seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11_quantified/counts_plot_table_noDoped.csv")

# list of round IDs
round_list <- unique(seqs_input$run_name)

# convert list of sequences into a matrix
seqs_matrix <- do.call(rbind, type.convert(strsplit(seqs_input$sequence, ""), as.is = TRUE))

# trim the sequence to keep the bases of the overhang
overhang_matrix <- seqs_matrix[,16:25]

# set overhang data length
data_length <- nrow(seqs_input)

# initialize data frame for base counts
overhang_data <- data.frame(
  run_ID = rep(NA, data_length),
  overhang = rep(NA, data_length),
  identity = rep(NA, data_length)
)

# loop over each base
for (seq_num in 1:data_length) {
  # add run ID
  overhang_data$run_ID[seq_num] <- seqs_input$run_name[seq_num]
  # add overhang sequences
  overhang_data$overhang[seq_num] <- paste(overhang_matrix[seq_num,], collapse="")
  # count number of matched bases to the expected overhang
  num_match <- mapply(function(x, y) {
    len <- length(x)
    sum(x[1:len] == y[1:len])
  }, strsplit(overhang_data$overhang[seq_num], ''), strsplit(expected_overhang, ''))
  # add percent identity to expected overhang sequence
  overhang_data$identity[seq_num] <- 100*num_match/10
}
  
# vector of identity bins
identity_bins <- seq(0, 100, by = 10)

# initialize data frame for identity bin counts
overhang_counts <- data.frame(
  run_ID = rep(NA, 11),
  identity = rep(NA, 11),
  counts = rep(NA, 11),
  frac_abundance = rep(NA, 11),
  identity_color = rep(NA, 11)
)
overhang_counts_out <- data.frame()

# loop over each run
for (run_num in min(round_list):max(round_list)) {
  # loop over identity bins
  for (bin_index in 1:11) {
    # add run ID
    overhang_counts$run_ID[bin_index] <- run_num
    # add identity
    overhang_counts$identity[bin_index] <- identity_bins[bin_index]
    # add overhang counts
    overhang_counts$counts[bin_index] <- nrow(overhang_data[overhang_data$identity == identity_bins[bin_index] & overhang_data$run_ID == run_num,])
    # add fraction abundance
    overhang_counts$frac_abundance[bin_index] <- overhang_counts$counts[bin_index]/quality[run_num]
    # add cluster plotting color
    overhang_counts$identity_color[bin_index] <- safe_colors[bin_index]
  }
  # add current run data
  overhang_counts_out <- rbind(overhang_counts_out, overhang_counts)
}

# add percent counts
overhang_counts_out$perc_abundance <- 100*overhang_counts_out$frac_abundance

# change zeros to NAs for plotting
overhang_counts_out$counts_na <- ifelse(overhang_counts_out$counts == 0, NA, overhang_counts_out$counts)
overhang_counts_out$perc_abundance_na <- ifelse(overhang_counts_out$perc_abundance == 0, NA, overhang_counts_out$perc_abundance)
overhang_counts_out$frac_abundance_na <- ifelse(overhang_counts_out$frac_abundance == 0, NA, overhang_counts_out$frac_abundance)

# subset the overhang counts by identity
#overhang_counts_subset <- overhang_counts_out[overhang_counts_out$identity >= 80,]

# create heatmap of overhang identity log counts
base_counts_plot <- ggplot(data = overhang_counts_out, aes(reorder(as.character(run_ID), run_ID), reorder(as.character(identity), identity), fill = log(counts_na))) + 
  theme_bw() +
  geom_tile(colour = "black") +
  ylab("Identity") +
  xlab("Round Number") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(name = "Log Counts",
                       low = safe_colors[3],
                       mid = safe_colors[4],
                       high = safe_colors[5],
                       midpoint = max(log(overhang_counts_out$counts))/2,
                       na.value = "white")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_log_counts.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(base_counts_plot)
dev.off()

# create heatmap of overhang identity fraction abundances
base_counts_plot <- ggplot(data = overhang_counts_out, aes(reorder(as.character(run_ID), run_ID), reorder(as.character(identity), identity), fill = frac_abundance_na)) + 
  theme_bw() +
  geom_tile(colour = "black") +
  ylab("Identity") +
  xlab("Round Number") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(name = "FA",
                       low = safe_colors[3],
                       mid = safe_colors[4],
                       high = safe_colors[5],
                       midpoint = max(overhang_counts_out$frac_abundance)/2,
                       na.value = "white")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_fraction_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(base_counts_plot)
dev.off()

# create heatmap of overhang identity percent
base_counts_plot <- ggplot(data = overhang_counts_out, aes(reorder(as.character(run_ID), run_ID), reorder(as.character(identity), identity), fill = perc_abundance_na)) + 
  theme_bw() +
  geom_tile(colour = "black") +
  ylab("Identity") +
  xlab("Round Number") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(name = "PA",
                       low = safe_colors[3],
                       mid = safe_colors[4],
                       high = safe_colors[5],
                       midpoint = max(overhang_counts_out$perc_abundance)/2,
                       na.value = "white")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_percent_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(base_counts_plot)
dev.off()

# create line plot of overhang identity percent
base_counts_plot <- ggplot(data=overhang_counts_out, aes(x=as.character(run_ID), y=perc_abundance, group=identity, color=identity_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Identity", labels = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), breaks = safe_colors[1:10], guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_percent_abundance_lines.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(base_counts_plot)
dev.off()

# export plotting data
write.csv(overhang_counts_out, file = paste(out_dir, "/above9_overhang_conservation.csv", sep = ""), row.names = FALSE, quote = FALSE)
