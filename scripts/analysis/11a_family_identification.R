#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# library
library(stringr)

# turn of scientific notation
options(scipen=10000)

# set the sequence length
seqLength <- 40

# set the run name
#runName <- "r1_S1_L001"
runName <- args[1]

# read in cluster family sequence data
#r8_peaks <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized/r8_S8_L001_cluster_peaks_table.csv")
r8_peaks <- read.csv(args[2])

# read in sequence count data
#seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))
#seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/run2/09a_quantified/noDoped/r1_S1_L001_counts_plot_table.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))
seqs_input <- read.csv(args[3], colClasses=c("run_name"="character", "counts_run_name"="character"))

# remove NAs
seqs_input <- na.omit(seqs_input)

# clean up round names
#seqs_input$counts_run_name <- str_split_i(str_split_i(seqs_input$counts_run_name, "_", 1), "", 2)

# set outputs directory
out_dir <- args[4]

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# setup data frame length
#data_length <- nrow(seqs_input) * nrow(r8_peaks)

# data frames of sequence data
seqs_data <- data.frame(
  run_name = "NA",
  sequence_ID = "NA", 
  sequence = "NA", 
  peak_cluster_ID = "NA", 
  peak_identity = "NA"
)
seqs_data_out <- data.frame()

# testing
#seqs_input <- seqs_input[seqs_input$run_name == 2  & seqs_input$counts_run_name == 2,]
#nrow(seqs_input)
#nrow(seqs_input[unique(seqs_input$sequence_ID),])

# loop over each sequence and compare with the peak to note if it has >= 90% similarity
for (seq_num in 1:nrow(seqs_input)) {
  # loop over each peak sequence
  for (peak_num in 1:nrow(r8_peaks)) {
    # compare the current sequence with the current peak
    numMatch <- mapply(
      function(x, y) {
        len <- length(x)
        sum(x[1:len] == y[1:len])
      }, 
      strsplit(seqs_input$sequence[seq_num], ''), 
      strsplit(r8_peaks$sequence[peak_num], '')
    )
    # determine identity percent
    identity_perc <- 100*numMatch/seqLength
    # set the index
    #index <- ((seq_num-1)*nrow(r8_peaks)) + peak_num
    # record data
    #print(seqs_data)
    seqs_data$run_name[1] <- seqs_input[seqs_input$sequence_ID == seq_num, "run_name"]
    seqs_data$sequence_ID[1] <- seq_num
    seqs_data$sequence[1] <- seqs_input[seqs_input$sequence_ID == seq_num, "sequence"]
    seqs_data$peak_cluster_ID[1] <- peak_num
    seqs_data$peak_identity[1] <- identity_perc
    # add the latest data
    seqs_data_out <- rbind(seqs_data_out, seqs_data)
  }
}

# export data
write.csv(seqs_data_out, file = paste(out_dir, "/", runName, "_family_identities.csv", sep = ""), row.names = FALSE, quote = FALSE)

# check how many sequences have >= 90% identity to each peak
#for (cluster_num in 0:(nrow(r8_peaks)-1)) {
#  print(cluster_num)
#  print(nrow(seqs_data_out[seqs_data_out[seqs_data_out$peak_cluster_ID == cluster_num,]$peak_identity >= 90,]))
#}

# keep sequences that have >= 90% identity to any peak
seqs_out <- seqs_data_out[seqs_data_out$peak_identity >= 90,]

# export data
write.csv(seqs_out, file = paste(out_dir, "/", runName, "_family_identities_atLeast90.csv", sep = ""), row.names = FALSE, quote = FALSE)
