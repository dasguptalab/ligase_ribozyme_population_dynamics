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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/SF1_top10_overhang_conservation_unique"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), palette.colors(palette = "Okabe-Ito"))

# A pairs with U and G pairs with C, but G also pairs with U
# G (C, U)
# A (U)
# U (A, G)
# C (G)

# overhang sequence: AUUCCGCA
# overhang wobble positions: 2, 3, and 6
# complement wobble positions: 3, 6, and 7
# 3
# 3, 6
# 3, 7
# 3, 6, 7
# 6
# 6, 7
# 7

# set the overhang sequence
overhang <- rev(c("ATTCCGCA"))

# store the overhang in an array
overhang <- unlist(strsplit(overhang, ""))

# reverse the overhang for comparing
rev_overhang <- rev(overhang)

# set expected overhang (ATTCCGCA) reverse complement (TGCGGAAT) sequence
#complement_seq <- c("T","G","C","G","G","A","A","T","G","C")
#complement_seq <- c("TGCGGAATGC")
#complement_seq <- c("T","G","C","G","G")
complement_seq <- c("TGCGGAAT")

# store the expected complement in an array
complement_seq <- unlist(strsplit(complement_seq, ""))

# wobble positions (TG[C,T]GG[A,G][A,G]T)
# TG[T]GGAAT
# TG[T]GG[G]AT
# TG[T]GGA[G]T
# TG[T]GG[G][G]T
# TGCGG[G]AT
# TGCGG[G][G]T
# TGCGGA[G]T

# read in sequence data
seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09b_quantified/counts_plot_table_noDoped.csv")

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
  gap = rep(NA, seq_data_length),
  wobble = rep(NA, seq_data_length)
)

# initialize loop variable
sequence_ID_last <- "0_0"
wobble_flag <- NA

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
    # set the wobble flag
    complement_data$wobble[seq_num] <- complement_data$wobble[seq_num-1]
    # jump to the end of the loop and stop parsing the current sequence
    next
  }
  # set last sequence_ID
  sequence_ID_last <- seqs_input$sequence_ID[seq_num]
  # loop over each base of the sequence
  for (base_index in 1:sliding_length) {
    # reset window variables
    expected_complement <- complement_seq
    # set end index
    end_index <- base_index + (complement_length-1)
    # get the next 7 bases to create the 8bp sliding window
    slide_window <- seq_matrix[seq_num,base_index:end_index]
    # count number of matched bases to the expected overhang complement
    num_match <- mapply(
      function(window, complement) {
        len <- length(window)
        sum(window[1:len] == complement[1:len])
      }, 
      slide_window, 
      expected_complement
    )
    # count number of matched bases (including wobbles) to the expected overhang complement
    num_match_wobble <- mapply(
      function(window, rev_over) {
        # initialize counter
        counter <- 0
        # retrieve sequence length
        len <- length(window)
        # loop over each non gaped base of the window
        for (i in 1:len) {
          # check if the current base is a match
          # accounting for wobble
          # A pairs with U and G pairs with C, but G also pairs with U
          if (rev_over[i] == "A" & window[i] == "T") { # A (U)
            counter <- counter +1
          }else if (rev_over[i] == "C" & window[i] == "G"){ # C (G)
            counter <- counter +1
          }else if (rev_over[i] == "G" & (window[i] == "C" || window[i] == "T")){ # G (C, U)
            counter <- counter +1
          }else if (rev_over[i] == "T" & (window[i] == "A" || window[i] == "G")){ # U (A, G)
            counter <- counter +1
          }
        }
        # return the match count
        counter
      }, 
      slide_window, 
      rev_overhang
    )
    # determine percent identity to expected overhang complement
    window_identity <- 100*sum(num_match)/complement_length
    window_identity_wobble <- 100*sum(num_match_wobble)/complement_length
    # check if the wobble identity is higher
    if(window_identity_wobble > window_identity){
      # flag that the wobble is detected
      wobble_flag <- "yes"
      # set window identity to the higher wobble identity
      window_identity <- window_identity_wobble
      # else, wobble is not necessary
    }else{
      # flag that the wobble is not detected
      wobble_flag <- "no"
    }
    # check if window identity is not higher than last best
    if (window_identity < complement_data$identity[seq_num]) {
      # jump to the end of the loop and stop parsing the current window
      next
    }
    # check if identity is = to 100*8/8 = 100 and wobble is not detected
    if(window_identity == 100 && wobble_flag == "no") {
      # store the current window sequence as the complement
      complement_data$complement[seq_num] <- paste(seq_matrix[seq_num,base_index:end_index], collapse="")
      # add percent identity to expected overhang complement
      complement_data$identity[seq_num] <- window_identity
      # set the longest subset window identity
      complement_data$identity_subset[seq_num] <-  window_identity
      # flag that the current window does not have a gap
      complement_data$gap[seq_num] <-  "no"
      # set the wobble flag
      complement_data$wobble[seq_num] <- wobble_flag
      # break loop and stop parsing the current sequence
      break
      # check for gaps
    }else{
      # initialize subset length variable and mismatch flag
      subset_length <- 0
      subset_longest <- 0
      mismatch_flag <- 0
      # loop over each non gaped base of the window
      for (window_index in 1:complement_length) {
        # check if the current base is a match
        #if (slide_window[window_index] == expected_complement[window_index]) {
        # A pairs with U and G pairs with C, but G also pairs with U
        if (rev_overhang[window_index] == "A" & slide_window[window_index] == "T") { # A (U)
          # increment subset length
          subset_length <- subset_length+1
        }else if (rev_overhang[window_index] == "C" & slide_window[window_index] == "G"){ # C (G)
          # increment subset length
          subset_length <- subset_length+1
        }else if (rev_overhang[window_index] == "G" & (slide_window[window_index] == "C" || slide_window[window_index] == "T")){ # G (C, U)
          # increment subset length
          subset_length <- subset_length+1
        }else if (rev_overhang[window_index] == "T" & (slide_window[window_index] == "A" || slide_window[window_index] == "G")){ # U (A, G)
          # increment subset length
          subset_length <- subset_length+1
        }else{ # mismatch
          # flag mismatch
          mismatch_flag <- 1
        }
        # check if mismatch
        if (mismatch_flag == 1){
          # check if the current subset length is longest
          if (subset_length > subset_longest) {
            subset_longest <- subset_length
          }
          # reset subset length
          subset_length <- 0
          # reset mismatch flag
          mismatch_flag <- 0
        }else{
          # check if the current subset length is longest
          if (subset_length > subset_longest) {
            subset_longest <- subset_length
          }
        }
      }
      # set longest window subset identity
      subset_identity <- 100*subset_longest/complement_length
      # check if the identity of the current non gaped subset matches the total
      if (subset_identity == window_identity & window_identity >= complement_data$identity[seq_num]) {
        # store the current window sequence as the complement
        complement_data$complement[seq_num] <- paste(seq_matrix[seq_num,base_index:end_index], collapse="")
        # add percent identity to expected overhang complement
        complement_data$identity[seq_num] <- window_identity
        # add subset percent identity to expected overhang complement
        complement_data$identity_subset[seq_num] <- subset_identity
        # flag that the current window does not have a gap
        complement_data$gap[seq_num] <-  "no"
        # set the wobble flag
        complement_data$wobble[seq_num] <- wobble_flag
      }else if (subset_identity != window_identity & window_identity > complement_data$identity[seq_num]){
        # store the current window sequence as the complement
        complement_data$complement[seq_num] <- paste(seq_matrix[seq_num,base_index:end_index], collapse="")
        # add percent identity to expected overhang complement
        complement_data$identity[seq_num] <- window_identity
        # add subset percent identity to expected overhang complement
        complement_data$identity_subset[seq_num] <- subset_identity
        # flag that the current window does have a gap
        complement_data$gap[seq_num] <-  "yes"
        # set the wobble flag
        complement_data$wobble[seq_num] <- wobble_flag
      }
    }
  }
}

# vectors of bins (total, non gaped, gaped)
identity_bins <- unique(complement_data$identity)
#plot_bins <- c(paste(identity_bins, "T", sep = "_"), paste(identity_bins, "C", sep = "_"), paste(identity_bins, "G", sep = "_"), paste(identity_bins, "W", sep = "_"), paste(identity_bins, "WG", sep = "_"))
plot_bins <- c(paste(identity_bins, "T", sep = "_"), paste(identity_bins, "C", sep = "_"), paste(identity_bins, "G", sep = "_"))

# set data length
data_length <- length(plot_bins)

# keep sequences with at least a 3bp non gaped match (100*3/8 = 37.5)
#complement_data_subset <- complement_data[complement_data$identity_subset >= 37.5,]
complement_data_subset <- complement_data

# initialize data frame for identity bin counts
complement_counts <- data.frame(
  run_name = rep(NA, data_length),
  identity = rep(NA, data_length),
  type = rep(NA, data_length),
  counts_unique = rep(NA, data_length),
  frac_abundance_unique = rep(NA, data_length),
  identity_type_color = rep(NA, data_length),
  identity_color = rep(NA, data_length),
  identity_label = rep(NA, data_length)
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
    if (cur_type == "T") { # total (including wobble)
      # add overhang complement unique counts
      complement_counts$counts_unique[bin_index] <- nrow(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num,])
    }else if (cur_type == "C") { # non gaped (including wobble)
      # add overhang complement unique counts
      complement_counts$counts_unique[bin_index] <- nrow(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "no" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num,])
    }else if (cur_type == "G") { # gaped (including wobble)
      # add overhang complement unique counts
      complement_counts$counts_unique[bin_index] <- nrow(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "yes" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num,])
    }
    # add fraction abundance
    complement_counts$frac_abundance_unique[bin_index] <- complement_counts$counts_unique[bin_index]/nrow(complement_data_subset[complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == run_num,])
    # add cluster plotting color
    complement_counts$identity_type_color[bin_index] <- safe_colors[bin_index]
  }
  # add current run data
  complement_counts_out <- rbind(complement_counts_out, complement_counts)
}

# add percent counts
complement_counts_out$perc_abundance_unique <- 100*complement_counts_out$frac_abundance_unique

# get lists of identities and types
identity_list <- unique(complement_counts_out$identity)
type_list <- unique(complement_counts_out$type)

# sort the identity list
#identity_list <- as.character(sort(as.numeric(identity_list), decreasing = FALSE))

# list of identity labels
identity_labels <- c(7, 6, 8, 5)

# set identity colors for plotting
complement_counts_out[complement_counts_out$identity == identity_list, "identity_color"] <- safe_colors[1:length(identity_list)]

# set identity labels for plotting
complement_counts_out[complement_counts_out$identity == identity_list, "identity_label"] <- identity_labels

# sort the data for plotting
complement_counts_sorted <- complement_counts_out[order(complement_counts_out$identity_label, decreasing = TRUE),]

# subset counts by type
complement_counts_total <- complement_counts_sorted[complement_counts_sorted$type == "T",]
complement_counts_consecutive <- complement_counts_sorted[complement_counts_sorted$type == "C",]
complement_counts_gap <- complement_counts_sorted[complement_counts_sorted$type == "G",]

# create line plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_total, aes(fill=identity, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("#DDCC77", "#117733", "#CC6677", "#88CCEE"), labels = c(8, 5, 6, 7)) +
  labs(fill = "Matched Bases") +
  ylab("Proportion") +
  xlab("Round")
#base_counts_plot <- ggplot(data=complement_counts_total, aes(x=as.character(run_name), y=perc_abundance_unique, group=identity, color=identity_color))+
#  geom_line(size = 1) +
#  geom_point() +
#  theme_classic(base_size = 14) +
#  scale_color_identity(name = "Matched Bases", labels = complement_counts_total$identity_label, breaks = complement_counts_total$identity_color, guide = "legend") +
#  ylab("Proportion") +
#  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of non gaped overhang identity percent
base_counts_plot <- ggplot(complement_counts_consecutive, aes(fill=identity, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("#DDCC77", "#117733", "#CC6677", "#88CCEE"), labels = c(8, 5, 6, 7)) +
  labs(fill = "Matched Bases") +
  ylab("Proportion") +
  xlab("Round")
#base_counts_plot <- ggplot(data=complement_counts_consecutive, aes(x=as.character(run_name), y=perc_abundance_unique, group=identity, color=identity_color))+
#  geom_line(size = 1) +
#  geom_point() +
#  theme_classic(base_size = 14) +
#  scale_color_identity(name = "Matched Bases", labels = complement_counts_consecutive$identity_label, breaks = complement_counts_consecutive$identity_color, guide = "legend") +
#  ylab("Proportion") +
#  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_non_gaped.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create line plot of gaped overhang identity percent
#complement_counts_gap_subset <- complement_counts_gap[complement_counts_gap$identity == "75" | complement_counts_gap$identity == "62.5",]
base_counts_plot <- ggplot(complement_counts_gap, aes(fill=identity, y=perc_abundance_unique, x=as.character(run_name))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 14) +
  #scale_fill_manual(values = c("#CC6677", "#88CCEE"), labels = c(6, 7)) +
  scale_fill_manual(values = c("#DDCC77", "#117733", "#CC6677", "#88CCEE"), labels = c(8, 5, 6, 7)) +
  labs(fill = "Matched Bases") +
  ylab("Proportion") +
  xlab("Round")
#base_counts_plot <- ggplot(data=complement_counts_gap, aes(x=as.character(run_name), y=perc_abundance_unique, group=identity, color=identity_color))+
#  geom_line(size = 1) +
#  geom_point() +
#  theme_classic(base_size = 14) +
#  scale_color_identity(name = "Matched Bases", labels = complement_counts_gap$identity_label, breaks = complement_counts_gap$identity_color, guide = "legend") +
#  ylab("Proportion") +
#  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_gaped.png", sep = "")
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
write.csv(complement_data, file = paste(out_dir, "/overhang_data_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_counts_sorted, file = paste(out_dir, "/overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_data_similar, file = paste(out_dir, "/overhang_data_similar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_data_disimilar, file = paste(out_dir, "/overhang_data_disimilar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_data_disimilar_r8, file = paste(out_dir, "/overhang_data_round8_disimilar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
#write.csv(complement_data_similar_seqs, file = paste(out_dir, "/overhang_data_similar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
#write.csv(complement_data_disimilar_seqs, file = paste(out_dir, "/overhang_data_disimilar_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
#write.csv(complement_data_similar_seqs_unique, file = paste(out_dir, "/overhang_data_similar_unique_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
#write.csv(complement_data_disimilar_seqs_unique, file = paste(out_dir, "/overhang_data_disimilar_unique_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
