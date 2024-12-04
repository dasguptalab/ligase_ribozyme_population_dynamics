#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
#library(plyr)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/02_family_abundances"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in cluster family sequence data
r8_seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_families/r8_S8_L001_counts_plot_table.csv")

# read in sequence count data
seqs_counts <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11_quantified_top10/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))

# reverse complement the sequences
#seqs_counts$sequence <- rev(chartr("ATGC","TACG",seqs_counts$sequence))

# add the fraction abundances
seqs_counts$frac_abundance <- NA
for (run_num in 1:8) {
  seqs_counts[seqs_counts$counts_run_name == run_num, "frac_abundance"] <- seqs_counts[seqs_counts$counts_run_name == run_num, "counts"]/quality[run_num]
}

# add log values
seqs_counts$log_counts <- log(seqs_counts$counts)
seqs_counts$log_frac_abundance <- log(seqs_counts$frac_abundance)

# set infinite and NA values equal to zero
is.na(seqs_counts)<-sapply(seqs_counts, is.infinite)

# re-order the data for plotting
seqs_counts <- seqs_counts[order(seqs_counts$log_counts, decreasing=TRUE),]

# setup midpoint values for plotting
seqs_counts_noNA <- seqs_counts
seqs_counts_noNA[is.na(seqs_counts_noNA)] <- 0
mid_log_counts <- max(seqs_counts_noNA$log_counts)/2
mid_log_frac_abundance <- log(min(seqs_counts_noNA[seqs_counts_noNA$frac_abundance != 0,"frac_abundance"]))/2

# list of cluster IDs
cluster_list <- unique(r8_seqs_family$cluster_ID)
cluster_list_out <- cluster_list+1

# setup data frame length
data_length <- (max(cluster_list)+1)*8

# data frame of cluster abundances
cluster_abundances <- data.frame(
  cluster_ID = rep(NA, data_length), 
  run_name = rep(NA, data_length),
  cluster_color = rep(NA, data_length),
  counts = rep(NA, data_length),
  frac_abundance = rep(NA, data_length)
  #frac_abundance_log = rep(NA, data_length)
)

# calculate fraction abundance per round
for (cluster_num in 0:max(cluster_list)) {
  # loop over each run
  for (run_num in 1:8) {
    # set the index
    index <- run_num+(cluster_num*8)
    # update family number for publishing
    cluster_out <- cluster_num+1
    # add cluster ID
    cluster_abundances$cluster_ID[index] <- paste("family", cluster_out, sep="_")
    # add run name
    cluster_abundances$run_name[index] <- run_num
    # add cluster plotting color
    cluster_abundances$cluster_color[index] <- safe_colors[cluster_out]
    # setup the column name
    #curr_col <- paste("r", run_num, "_counts", sep="")
    # add counts
    cluster_abundances$counts[index] <- sum(r8_seqs_family[r8_seqs_family$cluster_ID == cluster_num & r8_seqs_family$counts_run_name == run_num, "counts"])
    # add fraction abundance
    cluster_abundances$frac_abundance[index] <- cluster_abundances$counts[index]/quality[run_num]
  }
}

# add percentage abundances
cluster_abundances$perc_abundance <- 100*cluster_abundances$frac_abundance

# setup length of counts
counts_length <- 13*8

# line plot with counts per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=counts, group=cluster_ID, color=cluster_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Family", labels = cluster_list_out, breaks = safe_colors, guide = "legend") +
  ylab("Fraction Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/family_counts.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# line plot with log counts per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=log(counts), group=cluster_ID, color=cluster_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Family", labels = cluster_list_out, breaks = safe_colors, guide = "legend") +
  ylab("Log Fraction Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/family_log_counts.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# line plot with fraction abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=frac_abundance, group=cluster_ID, color=cluster_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Family", labels = cluster_list_out, breaks = safe_colors, guide = "legend") +
  ylab("Fraction Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/family_fraction_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# line plot with percent abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=perc_abundance, group=cluster_ID, color=cluster_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Family", labels = cluster_list_out, breaks = safe_colors, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/family_percent_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# export plotting data
write.csv(cluster_abundances, file = paste(out_dir, "/family_abundances.csv", sep = ""), row.names = FALSE, quote = FALSE)
