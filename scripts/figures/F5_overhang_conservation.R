#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# import libraries
library(ggplot2)
library(scales)
library(rcartocolor)
library(stringr)

# set the input round number
#round_num <- "8"
round_num <- args[1]
round_name <- paste("r", round_num, "_S", round_num, "_L001", sep = "")

# set outputs directory
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/F5_overhang_conservation"
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE)

# read in sequence data
#seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))
#seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified_all/r8_S8_L001_in_r8_S8_L001_counts_plot_table.csv"
seqsFile <- args[3]
seqs_input <- read.csv(seqsFile, colClasses=c("run_name"="character", "counts_run_name"="character"))


# create line plot of total overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_total, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=identity_color))+
  geom_line(size = 1) +
  geom_point() +
  theme_classic() +
  scale_color_identity(name = "Matched Bases", labels = complement_counts_total$identity_label, breaks = complement_counts_total$identity_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_percent_abundance_total.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of consecutive overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_consecutive, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=identity_color))+
  geom_line(size = 1) +
  geom_point() +
  theme_classic() +
  scale_color_identity(name = "Matched Bases", labels = complement_counts_consecutive$identity_label, breaks = complement_counts_consecutive$identity_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_percent_abundance_consecutive.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of gaped overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_gap, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=identity_color))+
  geom_line(size = 1) +
  geom_point() +
  theme_classic() +
  scale_color_identity(name = "Matched Bases", labels = complement_counts_gap$identity_label, breaks = complement_counts_gap$identity_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_percent_abundance_gaped.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# subset complement data by similarity
complement_data_similar <- complement_data_subset
complement_data_disimilar <- complement_data[complement_data$identity < 37.5 | complement_data$identity_subset < 37.5,]

# subset the round 8 data
complement_data_disimilar_r8 <- complement_data_disimilar[complement_data_disimilar$run_name == "8" & complement_data_disimilar$counts_run_name == "8",]

# keep sequence info
#complement_data_similar_seqs <- complement_data_similar[!duplicated(complement_data_similar$sequence_ID),"sequence"]
#complement_data_disimilar_seqs <- complement_data_disimilar[!duplicated(complement_data_disimilar$sequence_ID),"sequence"]

# keep unique sequences
#complement_data_similar_seqs_unique <- unique(complement_data_similar_seqs)
#complement_data_disimilar_seqs_unique <- unique(complement_data_disimilar_seqs)

# export data
write.csv(complement_data, file = paste(out_dir, "/", round_name, "_above9_overhang_data_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_counts_sorted, file = paste(out_dir, "/", round_name, "_above9_overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_data_similar, file = paste(out_dir, "/", round_name, "_above9_overhang_data_similar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_data_disimilar, file = paste(out_dir, "/", round_name, "_above9_overhang_data_disimilar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_data_disimilar_r8, file = paste(out_dir, "/", round_name, "_above9_overhang_data_round8_disimilar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
#write.csv(complement_data_similar_seqs, file = paste(out_dir, "/", round_name, "_above9_overhang_data_similar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
#write.csv(complement_data_disimilar_seqs, file = paste(out_dir, "/", round_name, "_above9_overhang_data_disimilar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
#write.csv(complement_data_similar_seqs_unique, file = paste(out_dir, "/", round_name, "_above9_overhang_data_similar_unique_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
#write.csv(complement_data_disimilar_seqs_unique, file = paste(out_dir, "/", round_name, "_above9_overhang_data_disimilar_unique_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
