#!/usr/bin/env Rscript

# R script to create overhang conservation plots
# usage: 05_F5_FS6_T3_TS3_TS4_overhang_conservation.R

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
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
#unique_reads <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)
#above_two_reads <- c(18, 19, 26, 27, 1585, 10626, 7230, 6315)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F5_overhang_conservation_above2"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F5_overhang_conservation_all"
#out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE)

# read in identity data
identityFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13b_overhang_conservation_above2/overhang_data_wobble.csv"
#identityFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13a_overhang_conservation_all/overhang_data_wobble.csv"
#seqsFile <- args[3]
complement_data_rounds <- read.csv(identityFile)

# replace NAs with zeros
complement_data_rounds[is.na(complement_data_rounds)] <- 0

# copy data frame for later plotting
complement_data <- complement_data_rounds

#keep sequences with at least a 4bp consecutive match (100*4/8 = 50)
#complement_data_similar <- complement_data[complement_data$identity >= 50, c("sequence", "counts", "counts_run_name")]
complement_data_dissimilar <- complement_data[complement_data$identity < 50, c("sequence", "counts", "counts_run_name")]
complement_data_dissimilar_round8 <- complement_data_dissimilar[complement_data_dissimilar$counts_run_name == "r8_S8_L001", c("sequence", "counts")]

# output dissimilar sequence data
write.csv(complement_data_dissimilar, file = paste(out_dir, "/overhang_conservation_wobble_dissimilar.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_data_dissimilar_round8, file = paste(out_dir, "/overhang_conservation_wobble_dissimilar_round8.csv", sep = ""), row.names = FALSE, quote = FALSE)

# vectors of bins (total, consecutive, gaped)
#tag_bins <- unique(complement_data$tag)
tag_bins <- c("0", "3", "4", "5", "3_3", "6", "7", "8")

# set data lengths
data_length <- length(tag_bins)

# list of round IDs
round_list <- unique(complement_data$run_name)

# initialize data frame for identity bin counts
complement_counts_sorted <- data.frame(
  run_name = rep(NA, data_length),
  tag = rep(NA, data_length),
  counts = rep(NA, data_length),
  counts_unique = rep(NA, data_length),
  frac_abundance = rep(NA, data_length),
  frac_abundance_unique = rep(NA, data_length)
)
complement_counts_out <- data.frame()

# loop over each run
for (run_num in min(round_list):max(round_list)) {
  # loop over tag bins
  for (bin_index in 1:data_length) {
    # retrieve current tag
    cur_tag <- tag_bins[bin_index]
    # add run ID
    complement_counts_sorted$run_name[bin_index] <- run_num
    # add tag
    complement_counts_sorted$tag[bin_index] <- cur_tag
    # set the counts run name
    counts_run_num <- paste("r", run_num, "_S", run_num, "_L001", sep = "")
    # set the number of unique sequences
    seq_data_length <- nrow(complement_data[complement_data$run_name == run_num,])
    # add overhang complement counts
    complement_counts_sorted$counts[bin_index] <- sum(complement_data[complement_data$tag == cur_tag & complement_data$run_name == run_num & complement_data$counts_run_name == counts_run_num, "counts"])
    # add overhang complement unique counts
    complement_counts_sorted$counts_unique[bin_index] <- nrow(complement_data[complement_data$tag == cur_tag & complement_data$run_name == run_num & complement_data$counts_run_name == counts_run_num,])
    # add fraction abundance
    complement_counts_sorted$frac_abundance_unique[bin_index] <- complement_counts_sorted$counts_unique[bin_index]/seq_data_length
    complement_counts_sorted$frac_abundance[bin_index] <- complement_counts_sorted$counts[bin_index]/quality[run_num]
    #complement_counts_sorted$frac_abundance[bin_index] <- complement_counts_sorted$counts[bin_index]/above_two_reads[run_num]
  }
  # add current run data
  complement_counts_out <- rbind(complement_counts_out, complement_counts_sorted)
}

# add percent counts
complement_counts_out$perc_abundance_unique <- 100*complement_counts_out$frac_abundance_unique
complement_counts_out$perc_abundance <- 100*complement_counts_out$frac_abundance

# export data
write.csv(complement_counts_out, file = paste(out_dir, "/", counts_run_num, "_overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)

# loop over each run
#for (run_num in 1:8) {
# calculate fraction abundance
#  complement_counts_sorted[complement_counts_sorted$run_name == run_num, "frac_abundance_unique"] <- complement_counts_sorted[complement_counts_sorted$run_name == run_num, "counts_unique"]/unique_reads[run_num]
#}
# calculate percent abundance
#complement_counts_sorted$perc_abundance_unique <- complement_counts_sorted$frac_abundance_unique*100

# get lists of identities and types
#identity_list <- unique(complement_counts_out$identity)
#identity_list <- sort(as.numeric(identity_list), decreasing=TRUE)
identity_list <- c(100, 87.5, 75, 75, 62.5, 50, 37.5, 0)
#tag_list <- unique(complement_counts_out$tag)
tag_list <- c("8", "7", "6", "3_3", "5", "4", "3", "0")

# add mapping table
identity_mappings <- data.frame(
  tag = tag_list,
  identity = identity_list,
  colors = safe_colors[1:8]
)

# add placeholder columns
complement_counts_out$colors <- NA

# loop over each possible tag label
for (label_num in 1:nrow(identity_mappings)) {
  # set tag colors for plotting
  ifelse(
    complement_counts_out$tag == identity_mappings$tag[label_num], 
    complement_counts_out[complement_counts_out$tag == identity_mappings$tag[label_num], "colors"] <- identity_mappings$colors[label_num], 
    NA
  )
}

# setup data for plotting
complement_counts_out$tag <- as.factor(complement_counts_out$tag)
#complement_counts_out <- complement_counts_out[complement_counts_out$tag != 3,]
#complement_counts_out <- complement_counts_out[complement_counts_out$tag != 4,]

# sort the data for plotting
complement_counts_total <- complement_counts_out[order(complement_counts_out$tag, decreasing = TRUE),]

## plots using sequencing read counts

# create line plot of total overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_total, aes(x=as.character(run_name), y=perc_abundance, group=tag, color=colors))+
  geom_line(size = 1.25) +
  geom_point(size = 2.25) +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Complementary\nBases", labels = complement_counts_total$tag, breaks = complement_counts_total$colors, guide = "legend") +
  scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10)) +
  guides(y = guide_axis(cap = "upper"), x = guide_axis(cap = "upper")) +
  ylab("Abundance") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_total.png", sep = "")
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

# view the exact complements
#View(complement_data[grep("100", complement_data$all_identities), ])

# import t0 data
t0File <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_t0_run1/overhang_conservation_t0.csv"
#t0File <- args[4]
t0_complement <- read.csv(t0File)
colnames(t0_complement) <- c("tag", "counts_unique", "colors", "perc_abundance_unique")
t0_complement$run_name <- "0"

# subset to necessary columns
complement_counts <- complement_counts_total[c("tag", "run_name", "perc_abundance_unique", "colors")]
complement_counts_t0 <- t0_complement[c("tag", "run_name", "perc_abundance_unique", "colors")]

# replace NAs with 0s
complement_counts[is.na(complement_counts)] <- 0

# remove factors
#complement_counts$bases <- as.numeric(as.character(complement_counts$bases))

# add t0 data
complement_counts_combined <- rbind(complement_counts, complement_counts_t0)

# add factors for plotting
#complement_counts_combined$bases <- as.factor(complement_counts_combined$bases)

# create bar plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_combined, aes(fill=tag, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_combined$tag), values = unique(complement_counts_combined$colors), labels = unique(complement_counts_combined$tag)) +
  #scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper")) +#, x = guide_axis(cap = "upper")) +
  labs(fill = "Complementary\nBases") +
  ylab("Proportion") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total_chart_t0.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# read in r0 overhang data
identityFile_t0 <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/overhang_conservation_t0_run1/r0_S0_L001_overhang_data_wobble.csv"
complement_data_t0 <- read.csv(identityFile_t0)

# replace NAs with zeros
complement_data_t0[is.na(complement_data_t0)] <- 0

# update t0 run name
complement_data_t0$run_name <- 0

# combine overhang data
complement_data_combined <- rbind(complement_data_t0[c("run_name","all_identities")], complement_data[c("run_name","all_identities")])

# initialize data column
complement_data_combined$num_regions <- NA
complement_data_combined$num_exact <- NA

# loop over each sequence
#for (seq_num in 1:nrow(complement_data_combined)) {
  # count number of complementary regions >= 37.5
  #complement_data_combined$num_regions[seq_num] <- str_count(complement_data_combined$all_identities[seq_num],";")+1
  # count number of exact matching regions
  #complement_data_combined$num_exact[seq_num] <- str_count(complement_data_combined$all_identities[seq_num],"100")+1
#}

# count number of complementary regions >= 37.5
complement_data_combined$num_regions <- str_count(complement_data_combined$all_identities,";")+1
# count number of exact matching regions
complement_data_combined$num_exact <- str_count(complement_data_combined$all_identities,"100")+1
# update for no matches
complement_data_combined[complement_data_combined$all_identities == 0, "num_regions"] <- 0
complement_data_combined[complement_data_combined$all_identities == 0, "num_exact"] <- 0

# get min and max number of regions 
min_regions <- min(complement_data_combined$num_regions, na.rm = TRUE)
max_regions <- max(complement_data_combined$num_regions, na.rm = TRUE)

# setup data lengths
data_length <- (max_regions+1)*(num_rounds+1)
num_t0 <- 1500000
min_rounds <- 0

# determine the frequency of the numbers of regions
region_data <- data.frame(
  run_name = rep(seq(0, num_rounds), (max_regions+1)),
  num_regions = c(rep(seq(0, max_regions), each = (num_rounds+1))),
  freq_regions = rep(NA, data_length),
  prop_regions = rep(NA, data_length),
  region_color = rep(NA, data_length) 
)

# loop over each round
for (round_num in min_rounds:num_rounds) {
  # loop over each number of regions
  for (region_num in min_regions:max_regions) {
    # set the freq of regions
    region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "freq_regions"] <- unname(table(complement_data_combined[complement_data_combined$run_name == round_num, "num_regions"]))[region_num+1]
    # check the round number
    if (round_num == 0) {
      # determine the proportion of sequences
      region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "prop_regions"] <- region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "freq_regions"]/num_t0
    }else{
      # determine the proportion of sequences
      region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "prop_regions"] <- region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "freq_regions"]/nrow(complement_data_combined[complement_data_combined$run_name == round_num,])
    }
    # set the region color
    region_data[region_data$run_name == round_num & region_data$num_regions == region_num, "region_color"] <- c(safe_colors[6], "#D55E00", "#0072B2", "#CC79A7", "#F0E442", "#009E73", "#56B4E9")[region_num+1]
  }
}
region_data$prop_regions <- region_data$prop_regions*100

# set NAs to 0
region_data[is.na(region_data)] <- 0

# loop over each round
freq_regions_rounds <- rep(NA, (num_rounds+1))
prop_regions_rounds <- rep(NA, (num_rounds+1))
for (round_num in 0:num_rounds) {
  # combine rows for number of regions >5
  freq_regions_rounds[round_num+1] <- sum(region_data[region_data$run_name == round_num & region_data$num_regions >= 6, "freq_regions"])
  prop_regions_rounds[round_num+1] <- sum(region_data[region_data$run_name == round_num & region_data$num_regions >= 6, "prop_regions"])
  
}

# create updated data frame for plotting
region_data_out <- region_data[region_data$num_regions <= 5,]
region_data_subset <- cbind(
  run_name = region_data[region_data$num_regions == 6, "run_name"],
  num_regions = rep(">5", (num_rounds+1)),
  freq_regions = freq_regions_rounds,
  prop_regions = prop_regions_rounds,
  region_color = rep(region_data[region_data$num_regions == 6, "region_color"][1], (num_rounds+1))
)
region_data_out <- rbind(region_data_out, region_data_subset)

# specify order of bars (from top to bottom)
region_data_out$num_regions <- factor(region_data_out$num_regions, levels=c("0", "1", "2", "3", "4", "5", ">5"))

# create bar plot of total overhang identity percent
regions_counts_plot <- ggplot(region_data_out, aes(fill=num_regions, y=as.numeric(prop_regions), x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = region_data_out$num_regions, values = region_data_out$region_color, labels = region_data_out$num_regions) +
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
