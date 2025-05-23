#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# R script to identify concerved overhang sequences

# set the library paths for the CRC servers
#.libPaths("/afs/crc.nd.edu/user/e/ebrooks5/R/x86_64-pc-linux-gnu-library/4.4")

# turn of scientific notation
options(scipen=10000)

# import libraries
#library(ggplot2)
library(scales)
#library(rcartocolor)
library(stringr)

# set the input round number
#round_num <- "8"
round_num <- args[1]
round_name <- paste("r", round_num, "_S", round_num, "_L001", sep = "")

# set outputs directory
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/F5_overhang_conservation_all"
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE)

# read in sequence data
#seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))
#seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified_all/r8_S8_L001_in_r8_S8_L001_counts_plot_table.csv"
seqsFile <- args[3]
seqs_input <- read.csv(seqsFile, colClasses=c("run_name"="character", "counts_run_name"="character"))

# subset to the round data
seqs_input <- seqs_input[seqs_input$run_name == round_num & seqs_input$counts_run_name == round_name,]

# color blind safe plotting palette
#safe_colors <- c(carto_pal(name="Safe"), palette.colors(palette = "Okabe-Ito"))

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

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
#unique_reads <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)
#filtered_reads <- c(18, 19, 26, 27, 1585, 10626, 7230, 6315)

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
  wobble = rep(NA, seq_data_length),
  location = rep(NA, seq_data_length),
  all_locations = rep(NA, seq_data_length),
  all_identities = rep(NA, seq_data_length)
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
    # set the location
    complement_data$location[seq_num] <- complement_data$location[seq_num-1]
    # set all complementary locations
    complement_data$all_locations[seq_num] <- complement_data$all_locations[seq_num-1]
    # set all complementary identities
    complement_data$all_identities[seq_num] <- complement_data$all_identities[seq_num-1]
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
        # loop over each consecutive base of the window
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
      # set the location
      complement_data$location[seq_num] <- paste(base_index, end_index, sep = "-")
      # update all complementary locations
      complement_data$all_locations[seq_num] <- paste(complement_data$all_locations[seq_num], paste(base_index, end_index, sep = "-"), sep = ";")
      # update all complementary identities
      complement_data$all_identities[seq_num] <- paste(complement_data$all_identities[seq_num], window_identity, sep = ";")
      # break loop and stop parsing the current sequence
      break
      # check for gaps
    }else{
      # how is binding influenced by matching to the overhang?
      # what about the 3' to 5' orientation? <- doesn't matter... same either way
      # initialize subset length variable and mismatch flag
      subset_length <- 0
      subset_longest <- 0
      mismatch_flag <- 0
      # loop over each consecutive base of the window
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
        # set the wobble flag
        complement_data$wobble[seq_num] <- wobble_flag
        # set the location
        complement_data$location[seq_num] <- paste(base_index, end_index, sep = "-")
        # check if the subset identity is at least 3/8
        if (subset_identity >= 37.5) {
          # update all complementary locations
          complement_data$all_locations[seq_num] <- paste(complement_data$all_locations[seq_num], paste(base_index, end_index, sep = "-"), sep = ";")
          # update all complementary identities
          complement_data$all_identities[seq_num] <- paste(complement_data$all_identities[seq_num], window_identity, sep = ";")
        }
      #}else if (subset_identity != window_identity & window_identity > complement_data$identity[seq_num]){
      }else if (subset_identity < window_identity & window_identity > complement_data$identity[seq_num]){
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
        # check if the subset identity is at least 3/8
        if (subset_identity >= 37.5) {
          # update all complementary locations
          complement_data$all_locations[seq_num] <- paste(complement_data$all_locations[seq_num], paste(base_index, end_index, sep = "-"), sep = ";")
          # update all complementary identities
          complement_data$all_identities[seq_num] <- paste(complement_data$all_identities[seq_num], window_identity, sep = ";")
        }
      }
    }
  }
}

# remove initializing NAs
complement_data$all_locations <- gsub("NA;", "", complement_data$all_locations)
complement_data$all_identities <- gsub("NA;", "", complement_data$all_identities)

# vectors of bins (total, consecutive, gaped)
identity_bins <- unique(complement_data$identity)
#plot_bins <- c(paste(identity_bins, "T", sep = "_"), paste(identity_bins, "C", sep = "_"), paste(identity_bins, "G", sep = "_"), paste(identity_bins, "W", sep = "_"), paste(identity_bins, "WG", sep = "_"))
plot_bins <- c(paste(identity_bins, "T", sep = "_"), paste(identity_bins, "C", sep = "_"), paste(identity_bins, "G", sep = "_"))

# set data length
data_length <- length(plot_bins)

# keep sequences with at least a 3bp consecutive match (100*3/8 = 37.5)
complement_data_subset <- complement_data[complement_data$identity_subset >= 37.5,]
#complement_data_subset <- complement_data

# initialize data frame for identity bin counts
complement_counts <- data.frame(
  run_name = rep(NA, data_length),
  identity = rep(NA, data_length),
  type = rep(NA, data_length),
  counts = rep(NA, data_length),
  counts_unique = rep(NA, data_length),
  frac_abundance = rep(NA, data_length),
  frac_abundance_unique = rep(NA, data_length)
  #identity_type_color = rep(NA, data_length),
  #identity_color = rep(NA, data_length),
  #identity_label = rep(NA, data_length)
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
    # set the counts run name
    counts_run_num <- paste("r", run_num, "_S", run_num, "_L001", sep = "")
    # check the type of complement
    if (cur_type == "T") { # total (including wobble)
      # add overhang complement counts
      complement_counts$counts[bin_index] <- sum(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == counts_run_num, "counts"])
      # add overhang complement unique counts
      complement_counts$counts_unique[bin_index] <- nrow(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == counts_run_num,])
    }else if (cur_type == "C") { # consecutive (including wobble)
      # add overhang complement counts
      complement_counts$counts[bin_index] <- sum(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "no" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == counts_run_num, "counts"])
      # add overhang complement unique counts
      complement_counts$counts_unique[bin_index] <- nrow(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "no" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == counts_run_num,])
    }else if (cur_type == "G") { # gaped (including wobble)
      # add overhang complement counts
      complement_counts$counts[bin_index] <- sum(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "yes" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == counts_run_num, "counts"])
      # add overhang complement unique counts
      complement_counts$counts_unique[bin_index] <- nrow(complement_data_subset[complement_data_subset$identity == cur_identity & complement_data_subset$gap == "yes" & complement_data_subset$run_name == run_num & complement_data_subset$counts_run_name == counts_run_num,])
    }
    # add fraction abundance
    #complement_counts$frac_abundance_unique[bin_index] <- complement_counts$counts_unique[bin_index]/filtered_reads[run_num]
    complement_counts$frac_abundance_unique[bin_index] <- complement_counts$counts_unique[bin_index]/seq_data_length
    complement_counts$frac_abundance[bin_index] <- complement_counts$counts[bin_index]/quality[run_num]
    # add identity plotting color
    #complement_counts$identity_type_color[bin_index] <- safe_colors[bin_index]
  }
  # add current run data
  complement_counts_out <- rbind(complement_counts_out, complement_counts)
}

# add percent counts
complement_counts_out$perc_abundance_unique <- 100*complement_counts_out$frac_abundance_unique
complement_counts_out$perc_abundance <- 100*complement_counts_out$frac_abundance

# get lists of identities and types
identity_list <- unique(complement_counts_out$identity)
identity_list <- sort(as.numeric(identity_list), decreasing=TRUE)
type_list <- unique(complement_counts_out$type)

# add mapping table
identity_mappings <- data.frame(
  identity = c(100, 87.5, 75, 62.5, 50, 37.5),
  bases = c(8, 7, 6, 5, 4, 3)
)

# list of identity labels
identity_labels <- identity_mappings[identity_mappings$identity %in% identity_list, "bases"]

# set identity colors for plotting
#complement_counts_out[complement_counts_out$identity %in% identity_list, "identity_color"] <- safe_colors[1:length(identity_list)]

# set identity labels for plotting
complement_counts_out[complement_counts_out$identity %in% identity_list, "identity_label"] <- identity_labels

# sort the data for plotting
complement_counts_sorted <- complement_counts_out[order(complement_counts_out$identity_label, decreasing = TRUE),]

# export data
write.csv(complement_data, file = paste(out_dir, "/", round_name, "_overhang_data_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_counts_sorted, file = paste(out_dir, "/", round_name, "_overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
