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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/05_above9_overhang_conservation"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), palette.colors(palette = "Okabe-Ito"))

# set expected overhang (ATTCCGCA) complement (TGCGGAAT) sequence
#complement_seq <- c("T","G","C","G","G","A","A","T","G","C")
#complement_seq <- c("TGCGGAATGC")
#complement_seq <- c("T","G","C","G","G")
complement_seq <- c("TGCGGAAT")

# store the expected complement in an array
expected_complement <- unlist(strsplit(complement_seq, ""))

# set overhang length
frag_length <- 5

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
#unique_reads <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)
above_9_reads <- c(5, 3, 5, 4, 283, 4001, 2703, 2100)

# read in sequence data
seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11_quantified_above9/counts_plot_table_noDoped.csv")

# list of round IDs
round_list <- unique(seqs_input$run_name)

# convert list of sequences into a matrix
seq_matrix <- do.call(rbind, type.convert(strsplit(seqs_input$sequence, ""), as.is = TRUE))

# trim the sequence to keep the bases of the overhang
#seq_matrix <- seqs_input_matrix[,16:20]

# set seqeunce and overhang complement lengths
seq_length <- 40
complement_length <- 8

# set length for window sliding
sliding_length <- seq_length - (complement_length-1)

# set overhang data length
seq_data_length <- nrow(seqs_input)

# initialize data frame for base counts
complement_data <- data.frame(
  run_name = rep(NA, seq_data_length),
  sequence_ID = rep(NA, seq_data_length),
  sequence = rep(NA, seq_data_length),
  counts = rep(NA, seq_data_length),
  counts_run_name = rep(NA, seq_data_length),
  complement = rep(NA, seq_data_length),
  identity = rep(0, seq_data_length),
  identity_subset = rep(0, seq_data_length),
  gap = rep(NA, seq_data_length)
)

# initialize loop variable
sequence_ID_last <- "0_0"

# loop over each sequence
for (seq_num in 1:seq_data_length) {
  # add run ID
  complement_data$run_name[seq_num] <- seqs_input$run_name[seq_num]
  # add sequence ID
  complement_data$sequence_ID[seq_num] <- seqs_input$sequence_ID[seq_num]
  # add sequence
  complement_data$sequence[seq_num] <- paste(seq_matrix[seq_num,], collapse="")
  # add counts
  complement_data$counts[seq_num] <- seqs_input$counts[seq_num]
  # add counts run name
  complement_data$counts_run_name[seq_num] <- seqs_input$counts_run_name[seq_num]
  # check if the current sequence has already been analyzed
  # if so, set the current sequence data to be the same as the previous and stop parsing
  if (seqs_input$sequence_ID[seq_num] == sequence_ID_last) {
    # set the window sequence
    complement_data$complement[seq_num] <- complement_data$complement[seq_num-1]
    # set the percent identity
    complement_data$identity[seq_num] <- complement_data$identity[seq_num-1]
    # set the longest subset window identity
    complement_data$identity_subset[seq_num] <-  complement_data$identity_subset[seq_num-1]
    # set window gap flag
    complement_data$gap[seq_num] <-  complement_data$gap[seq_num-1]
    # jump to the end of the loop and stop parsing the current sequence
    next
  }
  # set last sequence_ID
  sequence_ID_last <- seqs_input$sequence_ID[seq_num]
  # loop over each base of the sequence
  for (base_index in 1:sliding_length) {
    # set end index
    end_index <- base_index + (complement_length-1)
    # get the next 7 bases to create the 8bp sliding window
    slide_window <- seq_matrix[seq_num,base_index:end_index]
    # count number of matched bases to the expected overhang complement
    num_match <- mapply(
      function(x, y) {
        len <- length(x)
        sum(x[1:len] == y[1:len])
      }, 
      slide_window, 
      expected_complement
    )
    # determine percent identity to expected overhang complement
    window_identity <- 100*sum(num_match)/complement_length
    # check if window identity is not higher than last best
    if (window_identity < complement_data$identity[seq_num]) {
      # jump to the end of the loop and stop parsing the current window
      next
    }
    # check if identity is = to 100*8/8 = 100
    if(window_identity == 100) {
      # store the current window sequence as the complement
      complement_data$complement[seq_num] <- paste(seq_matrix[seq_num,base_index:end_index], collapse="")
      # add percent identity to expected overhang complement
      complement_data$identity[seq_num] <- window_identity
      # set the longest subset window identity
      complement_data$identity_subset[seq_num] <-  window_identity
      # flag that the current window does not have a gap
      complement_data$gap[seq_num] <-  "no"
      # break loop and stop parsing the current sequence
      break
    # check for gaps
    }else{
      # initialize subset length variable
      subset_length <- 0
      subset_longest <- 0
      # loop over each consecutive base of the window
      for (window_index in 1:complement_length) {
        # check if the current base is a match
        if (slide_window[window_index] == expected_complement[window_index]) {
          # increment subset length
          subset_length <- subset_length+1
          # check if current subset length is longest
          if (subset_length > subset_longest) {
            subset_longest <- subset_length
          }
        }else{ # mismatch
          # reset subset length
          subset_length <- 0
        }
      }
      # set longest window subset identity
      subset_identity <- 100*subset_longest/complement_length
      # check if the identity of the current consecutive subset matches the total
      if (subset_identity == window_identity & window_identity >= complement_data$identity[seq_num]) {
        # store the current window sequence as the complement
        complement_data$complement[seq_num] <- paste(seq_matrix[seq_num,base_index:end_index], collapse="")
        # add percent identity to expected overhang complement
        complement_data$identity[seq_num] <- window_identity
        # add subset percent identity to expected overhang complement
        complement_data$identity_subset[seq_num] <- subset_identity
        # flag that the current window does not have a gap
        complement_data$gap[seq_num] <-  "no"
      }else if (subset_identity != window_identity & window_identity > complement_data$identity[seq_num]){
        # store the current window sequence as the complement
        complement_data$complement[seq_num] <- paste(seq_matrix[seq_num,base_index:end_index], collapse="")
        # add percent identity to expected overhang complement
        complement_data$identity[seq_num] <- window_identity
        # add subset percent identity to expected overhang complement
        complement_data$identity_subset[seq_num] <- subset_identity
        # flag that the current window does have a gap
        complement_data$gap[seq_num] <-  "yes"
      }
    }
  }
}
  
# vectors of bins (total, consecutive, gaped)
identity_bins <- unique(complement_data$identity)
plot_bins <- c(paste(identity_bins, "T", sep = "_"), paste(identity_bins, "C", sep = "_"), paste(identity_bins, "G", sep = "_"))
  
# set data length
data_length <- length(plot_bins)

# keep sequences with at least a 3bp consecutive match (100*3/8 = 37.5)
complement_data_subset <- complement_data[complement_data$identity_subset >= 37.5,]

# initialize data frame for identity bin counts
complement_counts <- data.frame(
  run_name = rep(NA, data_length),
  identity = rep(NA, data_length),
  type = rep(NA, data_length),
  counts = rep(NA, data_length),
  counts_unique = rep(NA, data_length),
  frac_abundance = rep(NA, data_length),
  frac_abundance_unique = rep(NA, data_length),
  identity_type_color = rep(NA, data_length),
  identity_color = rep(NA, data_length)
)
complement_counts_out <- data.frame()

# To-do: consider recording beyond best scoring
# loop over each run
for (run_num in min(round_list):max(round_list)) {
  # loop over identity bins
  for (bin_index in 1:data_length) {
    # retrieve current identity
    cur_identity <- strsplit(plot_bins[bin_index], split = "_")[[1]][1]
    # retrieve current type
    cur_type <- strsplit(plot_bins[bin_index], split = "_")[[1]][2]
    # add run ID
    complement_counts$run_name[bin_index] <- run_num
    # add identity
    complement_counts$identity[bin_index] <- cur_identity
    # add type
    complement_counts$type[bin_index] <- cur_type
    # check the type of complement
    if (cur_type == "T") { # total
      # add overhang complement counts
      complement_counts$counts[bin_index] <- sum(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num, "counts"])
      # add overhang complement unique counts
      complement_counts$counts_unique[bin_index] <- nrow(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num,])
    }else if (cur_type == "C") { # consecutive
      # add overhang complement counts
      complement_counts$counts[bin_index] <- sum(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "no" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num, "counts"])
      # add overhang complement unique counts
      complement_counts$counts_unique[bin_index] <- nrow(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "no" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num,])
    }else if (cur_type == "G") { # gaped
      # add overhang complement counts
      complement_counts$counts[bin_index] <- sum(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "yes" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num, "counts"])
      # add overhang complement unique counts
      complement_counts$counts_unique[bin_index] <- nrow(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "yes" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num,])
    }
    # add fraction abundance
    complement_counts$frac_abundance_unique[bin_index] <- complement_counts$counts_unique[bin_index]/above_9_reads[run_num]
    complement_counts$frac_abundance[bin_index] <- complement_counts$counts[bin_index]/quality[run_num]
    # add cluster plotting color
    complement_counts$identity_type_color[bin_index] <- safe_colors[bin_index]
  }
  # add current run data
  complement_counts_out <- rbind(complement_counts_out, complement_counts)
}

# add percent counts
complement_counts_out$perc_abundance_unique <- 100*complement_counts_out$frac_abundance_unique
complement_counts_out$perc_abundance <- 100*complement_counts_out$frac_abundance

# get lists of identities and types
identity_list <- unique(complement_counts_out$identity)
type_list <- unique(complement_counts_out$type)

# set identity colors for plotting
complement_counts_out[complement_counts_out$identity == identity_list, "identity_color"] <- safe_colors[1:length(identity_list)]

# subset counts by type
complement_counts_total <- complement_counts_out[complement_counts_out$type == "T",]
complement_counts_consecutive <- complement_counts_out[complement_counts_out$type == "C",]
complement_counts_gap <- complement_counts_out[complement_counts_out$type == "G",]

# create line plot of total overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_total, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=identity_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Identity", labels = complement_counts_total$identity, breaks = complement_counts_total$identity_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_percent_abundance_reads_total.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(base_counts_plot)
dev.off()

# create line plot of consecutive overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_consecutive, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=identity_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Identity", labels = complement_counts_consecutive$identity, breaks = complement_counts_consecutive$identity_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_percent_abundance_reads_consecutive.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(base_counts_plot)
dev.off()

# create line plot of gaped overhang identity percent
base_counts_plot <- ggplot(data=complement_counts_gap, aes(x=as.character(run_name), y=perc_abundance, group=identity, color=identity_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Identity", labels = complement_counts_gap$identity, breaks = complement_counts_gap$identity_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/above9_overhang_percent_abundance_reads_gaped.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(base_counts_plot)
dev.off()

# export plotting data
write.csv(complement_counts_out, file = paste(out_dir, "/above9_overhang_conservation.csv", sep = ""), row.names = FALSE, quote = FALSE)
