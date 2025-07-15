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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/12_doped_peak_identity"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# set sequence lengths
seqLength <- 40

# set must abundant peak sequence from clustering round 8
peakSeq <- "GAATGCTGCCAACCGTGCGGGCTAATTGGCAGACTGAGCT"

# read in doped sequence data
seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/doped21-r3_S12_L001_table.csv")

# reverse complement the sequences
seqs_family$sequence <- rev(chartr("ATGC","TACG",seqs_family$sequence))

# filter sequences to those with >= 90% identity to the peak sequence
seqs_family$peak_ID <- mapply(function(x,y) sum(x==y),strsplit(peakSeq,""),strsplit(seqs_family$sequence,""))
seqs_family$peak_ID <- 100*seqs_family$peak_ID/seqLength

# check number of sequences with at least 90% ID to the peak
nrow(seqs_family[seqs_family$peak_ID >= 90,])

# line plot with identities
identity_plot <- ggplot(data=seqs_family, aes(x=peak_ID))+
  geom_bar() +
  theme_bw() +
  #ylab("Frequency") +
  xlab("Percent ID")
# save the plot
exportFile <- paste(out_dir, "/doped_peak_identities.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(identity_plot)
dev.off()

# export plotting data
write.csv(seqs_family, file = paste(out_dir, "/doped_peak_identities.csv", sep = ""), row.names = FALSE, quote = FALSE)

