#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# import libraries
library(ggplot2)
library(rcartocolor)
library(stringr)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), palette.colors(palette = "Okabe-Ito"))

# set the number of rounds
num_rounds <- 1

# numbers of high quality reads
#quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
#unique_reads <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)
unique_reads <- c(1500000)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_t0_run1"
#out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE)

# read in complement data
seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13e_conservation_t0/overhang_data_wobble_run1.csv"
#seqsFile <- args[2]
complement_data <- read.csv(seqsFile)

#keep sequences with at least a 3bp consecutive match (100*3/8 = 37.5)
complement_data_dissimilar <- complement_data[complement_data$identity_subset < 37.5, c("sequence")]

# output dissimilar sequence data
write.csv(complement_data_dissimilar, file = paste(out_dir, "/overhang_conservation_wobble_dissimilar_t0.csv", sep = ""), row.names = FALSE, quote = FALSE)

# clean up columns
#complement_counts_sorted$identity_color <- NULL
#complement_counts_sorted$identity_type_color <- NULL
#complement_counts_sorted$identity_label <- NULL

# get lists of identities and types
identity_list <- unique(complement_data$identity)
identity_list <- sort(as.numeric(identity_list), decreasing=TRUE)

# add mapping table
identity_mappings <- data.frame(
  identity = c(100, 87.5, 75, 62.5, 50, 37.5),
  bases = c(8, 7, 6, 5, 4, 3),
  colors = safe_colors[1:6]
)

# add placeholder columns
complement_data$bases <- NA
complement_data$colors <- NA

# loop over each possible identity label
for (label_num in 1:nrow(identity_mappings)) {
  # set identity labels for plotting
  ifelse(
    complement_data$identity == identity_mappings$identity[label_num], 
    complement_data[complement_data$identity == identity_mappings$identity[label_num], "bases"] <- identity_mappings$bases[label_num], 
    NA
  )
  ifelse(
    complement_data$identity == identity_mappings$identity[label_num], 
    complement_data[complement_data$identity == identity_mappings$identity[label_num], "colors"] <- identity_mappings$colors[label_num], 
    NA
  )
}

# determine the frequency of the numbers of bases
max_bases <- max(complement_data$bases, na.rm = TRUE)
min_bases <- min(complement_data$bases, na.rm = TRUE)
base_names <- seq(min_bases, max_bases)
base_tbl <- table(complement_data$bases)
base_data <- data.frame(
  num_bases = as.character(base_names),
  freq_bases = c(unname(base_tbl)),
  #base_color = rev(safe_colors[1:length(base_tbl)])
  base_color = rev(c(safe_colors[1:5], safe_colors[12]))
)

# basic piechart of region frequencies
base_counts_plot <- ggplot(base_data, aes(x="", y=freq_bases, fill=num_bases)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = base_data$base_color, guide = guide_legend(reverse = TRUE)) +
  labs(fill = "Compementary\nNucleotides")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_total_t0.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# loop over each sequence
#for (seq_num in 1:nrow(complement_data)) {
  # count number of complementary regions >= 37.5
#  complement_data$num_regions[seq_num] <- str_count(complement_data$all_identities[seq_num],";")+1
  # count number of exact matching regions
#  complement_data$num_exact[seq_num] <- str_count(complement_data$all_identities[seq_num],"100")+1
#}

# determine the frequency of the numbers of regions
#max_regions <- max(complement_data$num_regions, na.rm = TRUE)
#region_data <- data.frame(
#  num_regions = as.character(seq(1, max_regions)),
#  freq_regions = c(unname(table(complement_data$num_regions))),
#region_color = safe_colors[1:max_regions]
#  region_color = c(safe_colors[1:5], safe_colors[12])
#)

# basic piechart of region frequencies
#ggplot(region_data, aes(x="", y=freq_regions, fill=num_regions)) +
#  geom_bar(stat="identity", width=1, color="white") +
#  coord_polar("y", start=0) +
#  theme_void() +
#  scale_fill_manual(values = region_data$region_color) +
#  labs(fill = "Number of\nRegions")
