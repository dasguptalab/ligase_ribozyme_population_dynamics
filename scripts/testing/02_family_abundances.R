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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/02_family_abundances"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in sequence count data
seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11_quantified_above9/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))

# read in sequences that have at least 90% identity to any peak
seqs_identities <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/00a_family_identification_above9/family_identities_above9_atLeast90.csv")

# subset to keep round 8 data
#r8_seqs_90 <- seqs_identities[seqs_identities$run_name == 8,]
  
# initialize data frame
seqs_families <- data.frame()
peak_cluster_IDs <- NULL

# loop over each sequences that has at least 90% identity to any peak
for (seq_num in 1:nrow(seqs_identities)) {
  # keep count data for sequences that have at least 90% identity to any peak
  seqs_90_data <- seqs_input[seqs_input$sequence_ID == seqs_identities$sequence_ID[seq_num],]
  # add counts
  seqs_families <- rbind(seqs_families, seqs_90_data)
  # add peak cluster ID to vector
  peak_cluster_IDs <- c(peak_cluster_IDs, rep(seqs_identities$peak_cluster_ID[seq_num], nrow(seqs_90_data)))
}

# add peak cluster IDs
seqs_families <- cbind(seqs_families, peak_cluster_IDs)

# list of cluster IDs
fam_list_out <- seq(1, 13)

# list of cluster IDs in order of abundance in round 8
r8_fams <- data.frame(
  fam_ID = fam_list_out,
  fam_color = safe_colors,
  cluster_ID = c(1, 3, 0, 2, 5, 8, 4, 7, 11, 6, 10, 12, 9)
)

# setup data frame length
data_length <- (max(fam_list_out))*8

# data frame of cluster abundances
cluster_abundances <- data.frame(
  cluster_ID = rep(NA, data_length), 
  fam_ID = rep(NA, data_length),
  fam_color = rep(NA, data_length),
  run_name = rep(NA, data_length),
  counts = rep(NA, data_length),
  frac_abundance = rep(NA, data_length)
)

# calculate fraction abundance per round
for (cluster_num in 0:max(fam_list_out-1)) {
  # loop over each run
  for (run_num in 1:8) {
    # set the index
    index <- run_num+(cluster_num*8)
    # update family number for publishing
    cluster_out <- cluster_num+1
    # add cluster ID
    cluster_abundances$cluster_ID[index] <- cluster_num
    # add family number
    cluster_abundances$fam_ID[index] <- r8_fams[r8_fams$cluster_ID == cluster_num, "fam_ID"]
    # add family color
    cluster_abundances$fam_color[index] <- r8_fams[r8_fams$cluster_ID == cluster_num, "fam_color"]
    # add run name
    cluster_abundances$run_name[index] <- run_num
    # add counts
    cluster_abundances$counts[index] <- sum(seqs_families[seqs_families$run_name == run_num & seqs_families$counts_run_name == run_num & seqs_families$peak_cluster_IDs == cluster_num, "counts"])
    # add fraction abundance
    cluster_abundances$frac_abundance[index] <- cluster_abundances$counts[index]/quality[run_num]
  }
}

# add percentage abundances
cluster_abundances$perc_abundance <- 100*cluster_abundances$frac_abundance

# line plot with percent abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_abundances, aes(x=as.character(run_name), y=perc_abundance, group=fam_ID, color=fam_color))+
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_identity(name = "Family", labels = r8_fams$fam_ID, breaks = r8_fams$fam_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/family_percent_abundances_atLeast90.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(cluster_abundances_plot)
dev.off()

# export plotting data
write.csv(cluster_abundances, file = paste(out_dir, "/family_abundances_atLeast90.csv", sep = ""), row.names = FALSE, quote = FALSE)
