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
num_rounds <- 8

# numbers of high quality reads
#quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
unique_reads <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F5_overhang_conservation_above2"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F5_overhang_conservation_all"
#out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE)

# read in sequence data
seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13b_overhang_conservation_above2/overhang_conservation_wobble.csv"
#seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13a_overhang_conservation_all/overhang_conservation_wobble.csv"
#seqsFile <- args[3]
complement_counts_sorted <- read.csv(seqsFile)

# loop over each run
#for (run_num in 1:8) {
  # calculate fraction abundance
#  complement_counts_sorted[complement_counts_sorted$run_name == run_num, "frac_abundance_unique"] <- complement_counts_sorted[complement_counts_sorted$run_name == run_num, "counts_unique"]/unique_reads[run_num]
#}
# calculate percent abundance
#complement_counts_sorted$perc_abundance_unique <- complement_counts_sorted$frac_abundance_unique*100

# read in identity data
identityFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13b_overhang_conservation_above2/overhang_data_wobble.csv"
#identityFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13a_overhang_conservation_all/overhang_data_wobble.csv"
#seqsFile <- args[3]
complement_data <- read.csv(identityFile)

#keep sequences with at least a 3bp consecutive match (100*3/8 = 37.5)
#complement_data_similar <- complement_data[complement_data$identity_subset >= 37.5, c("sequence", "counts", "counts_run_name")]
complement_data_dissimilar <- complement_data[complement_data$identity_subset < 37.5, c("sequence", "counts", "counts_run_name")]
complement_data_dissimilar_round8 <- complement_data_dissimilar[complement_data_dissimilar$counts_run_name == "r8_S8_L001", c("sequence", "counts")]

# output dissimilar sequence data
write.csv(complement_data_dissimilar, file = paste(out_dir, "/overhang_conservation_wobble_dissimilar.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_data_dissimilar_round8, file = paste(out_dir, "/overhang_conservation_wobble_dissimilar_round8.csv", sep = ""), row.names = FALSE, quote = FALSE)

# clean up columns
#complement_counts_sorted$identity_color <- NULL
#complement_counts_sorted$identity_type_color <- NULL
#complement_counts_sorted$identity_label <- NULL

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
#complement_counts_sorted <- complement_counts_sorted[complement_counts_sorted$bases != 3,]
#complement_counts_sorted <- complement_counts_sorted[complement_counts_sorted$bases != 4,]

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
  geom_point(size = 2.25) +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Complementary\nBases", labels = complement_counts_total$bases, breaks = complement_counts_total$colors, guide = "legend") +
  scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10), labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper"), x = guide_axis(cap = "upper")) +
  ylab("Abundance") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_total.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of consecutive overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_consecutive, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=colors))+
  geom_line(size = 1.25) +
  geom_point(size = 2.25) +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Complementary\nBases", labels = complement_counts_consecutive$bases, breaks = complement_counts_consecutive$colors, guide = "legend") +
  scale_y_continuous(limits=c(0, 30), breaks=seq(0, 30, 5), labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper"), x = guide_axis(cap = "upper")) +
  ylab("Abundance") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_non_gaped.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of gaped overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_gap, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=colors))+
  geom_line(size = 1.25) +
  geom_point(size = 2.25) +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Complementary\nBases", labels = complement_counts_gap$bases, breaks = complement_counts_gap$colors, guide = "legend") +
  scale_y_continuous(limits=c(0, 50), breaks=seq(0, 50, 10), labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper"), x = guide_axis(cap = "upper")) +
  ylab("Abundance") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_gaped.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create bar plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_total, aes(fill=bases, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_total$bases), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$bases)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper")) +#, x = guide_axis(cap = "upper")) +
  labs(fill = "Complementary\nBases") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total_chart.png", sep = "")
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

# initialize data column
complement_data$num_regions <- NA
complement_data$num_exact <- NA

# loop over each sequence
for (seq_num in 1:nrow(complement_data)) {
  # count number of complementary regions >= 37.5
  complement_data$num_regions[seq_num] <- str_count(complement_data$all_identities[seq_num],";")+1
  # count number of exact matching regions
  complement_data$num_exact[seq_num] <- str_count(complement_data$all_identities[seq_num],"100")+1
}

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

# get max number of regions 
max_regions <- max(complement_data$num_regions, na.rm = TRUE)

# setup data length
data_length <- max_regions*num_rounds

# determine the frequency of the numbers of regions
region_data <- data.frame(
  run_name = rep(seq(1, num_rounds), max_regions),
  num_regions = c(rep(1, num_rounds), rep(2, num_rounds), rep(3, num_rounds), rep(4, num_rounds), rep(5, num_rounds), rep(6, num_rounds)),
  freq_regions = rep(NA, data_length),
  prop_regions = rep(NA, data_length),
  region_color = rep(NA, data_length) 
)

# loop over each round
for (round_num in 1:num_rounds) {
  # loop over each number of regions
  for (region_num in 1:max_regions) {
    # set the freq of regions
    region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "freq_regions"] <- unname(table(complement_data[complement_data$run_name == round_num, "num_regions"]))[region_num]
    # determine the proportion of sequences
    region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "prop_regions"] <- region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "freq_regions"]/nrow(complement_data[complement_data$run_name == round_num,])
    # set the region color
    region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "region_color"] <- c(safe_colors[6], "#D55E00", "#0072B2", "#F0E442", "#009E73", "#999999")[region_num]
  }
}
region_data$prop_regions <- region_data$prop_regions*100
  
# set NAs to 0
region_data[is.na(region_data)] <- 0

# create bar plot of total overhang identity percent
regions_counts_plot <- ggplot(region_data, aes(fill=as.character(num_regions), y=prop_regions, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = region_data$num_regions, values = region_data$region_color, labels = region_data$num_regions) +
  #scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper")) +#, x = guide_axis(cap = "upper")) +
  labs(fill = "Number of\nRegions") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_regions_counts_total_chart.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(regions_counts_plot)
dev.off()

# view the exact complements
#View(complement_data[grep("100", complement_data$all_identities), ])

# import t0 data
t0File <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_t0_run1/overhang_conservation_t0.csv"
#t0File <- args[4]
t0_complement <- read.csv(t0File)
colnames(t0_complement) <- c("bases", "counts_unique", "colors", "perc_abundance_unique")
t0_complement$run_name <- "0"

# subset to necessary columns
complement_counts <- complement_counts_total[c("bases", "run_name", "perc_abundance_unique", "colors")]
complement_counts_t0 <- t0_complement[c("bases", "run_name", "perc_abundance_unique", "colors")]

# remove factors
complement_counts$bases <- as.numeric(as.character(complement_counts$bases))

# add t0 data
complement_counts_combined <- rbind(complement_counts, complement_counts_t0)

# add factors for plotting
complement_counts_combined$bases <- as.factor(complement_counts_combined$bases)

# create bar plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_combined, aes(fill=bases, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_combined$bases), values = unique(complement_counts_combined$colors), labels = unique(complement_counts_combined$bases)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper")) +#, x = guide_axis(cap = "upper")) +
  labs(fill = "Complementary\nBases") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total_chart_t0.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()
