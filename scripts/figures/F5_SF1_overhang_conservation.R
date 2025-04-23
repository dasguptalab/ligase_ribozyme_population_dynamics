#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# import libraries
library(ggplot2)
library(rcartocolor)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), palette.colors(palette = "Okabe-Ito"))

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/F5_SF1_overhang_conservation_above2"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/F5_SF1_overhang_conservation_all"
#out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE)

# read in sequence data
seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13b_overhang_conservation_above2/overhang_conservation_wobble.csv"
#seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13a_overhang_conservation_all/overhang_conservation_wobble.csv"
#seqsFile <- args[3]
complement_counts_sorted <- read.csv(seqsFile)

# clean up columns
complement_counts_sorted$identity_color <- NULL
complement_counts_sorted$identity_type_color <- NULL
complement_counts_sorted$identity_label <- NULL

# get lists of identities and types
identity_list <- unique(complement_counts_sorted$identity)
identity_list <- sort(as.numeric(identity_list), decreasing=TRUE)

# add mapping table
identity_mappings <- data.frame(
  identity = c(100, 87.5, 75, 62.5, 50, 37.5),
  bases = c(8, 7, 6, 5, 4, 3),
  colors = safe_colors[1:6]
)

# add placeholder columns
complement_counts_sorted$bases <- NA
complement_counts_sorted$colors <- NA

# loop over each possible identity label
for (label_num in 1:nrow(identity_mappings)) {
  # set identity labels for plotting
  ifelse(
    complement_counts_sorted$identity == identity_mappings$identity[label_num], 
    complement_counts_sorted[complement_counts_sorted$identity == identity_mappings$identity[label_num], "bases"] <- identity_mappings$bases[label_num], 
    NA
  )
  ifelse(
    complement_counts_sorted$identity == identity_mappings$identity[label_num], 
    complement_counts_sorted[complement_counts_sorted$identity == identity_mappings$identity[label_num], "colors"] <- identity_mappings$colors[label_num], 
    NA
  )
}

# setup data for plotting
complement_counts_sorted$bases <- as.factor(complement_counts_sorted$bases)
complement_counts_sorted <- complement_counts_sorted[complement_counts_sorted$bases != 4,]

# subset counts by type
complement_counts_total <- complement_counts_sorted[complement_counts_sorted$type == "T",]
complement_counts_consecutive <- complement_counts_sorted[complement_counts_sorted$type == "C",]
complement_counts_gap <- complement_counts_sorted[complement_counts_sorted$type == "G",]

# sort the data for plotting
complement_counts_total <- complement_counts_total[order(complement_counts_total$bases, decreasing = TRUE),]
complement_counts_consecutive <- complement_counts_consecutive[order(complement_counts_consecutive$bases, decreasing = TRUE),]
complement_counts_gap <- complement_counts_gap[order(complement_counts_gap$bases, decreasing = TRUE),]

## plots using sequencing read counts

# create line plot of total overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_total, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=colors))+
  geom_line(size = 1.25) +
  geom_point() +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Matched Bases", labels = complement_counts_total$bases, breaks = complement_counts_total$colors, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_total.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of consecutive overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_consecutive, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=colors))+
  geom_line(size = 1.25) +
  geom_point() +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Matched Bases", labels = complement_counts_consecutive$bases, breaks = complement_counts_consecutive$colors, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_non_gaped.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of gaped overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_gap, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=colors))+
  geom_line(size = 1.25) +
  geom_point() +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Matched Bases", labels = complement_counts_gap$bases, breaks = complement_counts_gap$colors, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_gaped.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_total, aes(fill=bases, y=perc_abundance, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_total$bases), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$bases)) +
  labs(fill = "Matched Bases") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_total_chart.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of non gaped overhang identity percent
base_counts_plot <- ggplot(complement_counts_consecutive, aes(fill=bases, y=perc_abundance, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 14) +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_total$bases), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$bases)) +
  labs(fill = "Matched Bases") +
  ylab("Proportion") +
  xlab("Round")
exportFile <- paste(out_dir, "/overhang_percent_abundance_non_gaped_chart.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of gaped overhang identity percent
base_counts_plot <- ggplot(complement_counts_gap, aes(fill=bases, y=perc_abundance, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_total$bases), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$bases)) +
  labs(fill = "Matched Bases") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_gaped_chart.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

## plots using unique sequence counts

# create line plot of total overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_total, aes(x=as.character(run_name), y=perc_abundance_unique, group=identity, color=colors))+
  geom_line(size = 1.25) +
  geom_point() +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Matched Bases", labels = complement_counts_total$bases, breaks = complement_counts_total$colors, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of consecutive overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_consecutive, aes(x=as.character(run_name), y=perc_abundance_unique, group=identity, color=colors))+
  geom_line(size = 1.25) +
  geom_point() +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Matched Bases", labels = complement_counts_consecutive$bases, breaks = complement_counts_consecutive$colors, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_non_gaped.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of gaped overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_gap, aes(x=as.character(run_name), y=perc_abundance_unique, group=identity, color=colors))+
  geom_line(size = 1.25) +
  geom_point() +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Matched Bases", labels = complement_counts_gap$bases, breaks = complement_counts_gap$colors, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_gaped.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_total, aes(fill=bases, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_total$bases), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$bases)) +
  labs(fill = "Matched Bases") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total_chart.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of non gaped overhang identity percent
base_counts_plot <- ggplot(complement_counts_consecutive, aes(fill=bases, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 14) +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_total$bases), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$bases)) +
  labs(fill = "Matched Bases") +
  ylab("Proportion") +
  xlab("Round")
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_non_gaped_chart.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of gaped overhang identity percent
base_counts_plot <- ggplot(complement_counts_gap, aes(fill=bases, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_total$bases), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$bases)) +
  labs(fill = "Matched Bases") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_gaped_chart.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# to-do: double check
# subset complement data by similarity
#complement_data_similar <- complement_counts_sorted
#complement_data_disimilar <- complement_counts_sorted[complement_counts_sorted$identity < 37.5 | complement_counts_sorted$identity_subset < 37.5,]

# subset the round 8 data
#complement_data_disimilar_r8 <- complement_data_disimilar[complement_data_disimilar$run_name == "8" & complement_data_disimilar$counts_run_name == "8",]

# keep sequence info
#complement_data_similar_seqs <- complement_data_similar[!duplicated(complement_data_similar$sequence_ID),"sequence"]
#complement_data_disimilar_seqs <- complement_data_disimilar[!duplicated(complement_data_disimilar$sequence_ID),"sequence"]

# keep unique sequences
#complement_data_similar_seqs_unique <- unique(complement_data_similar_seqs)
#complement_data_disimilar_seqs_unique <- unique(complement_data_disimilar_seqs)

# export data
write.csv(complement_counts_sorted, file = paste(out_dir, "/overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_counts_total, file = paste(out_dir, "/overhang_conservation_wobble_total.csv", sep = ""), row.names = FALSE, quote = FALSE)
