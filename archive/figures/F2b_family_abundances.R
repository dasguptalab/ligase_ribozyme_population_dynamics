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
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/F2b_family_abundances"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
#unique <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)

# read in family data
fam_data <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/ST2_family_table/r8_family_count_data.csv")

# read in sequences that have at least 90% identity to any peak
seqs_identities <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11b_family_identification_above2/family_identities_max_atLeast90.csv")

# list of cluster and run IDs
cluster_list <- unique(seqs_identities$peak_cluster_ID)
run_list <- unique(seqs_identities$run_name)

# set the data length
data_length <- length(cluster_list)*length(run_list)
  
# create cluster data frame
cluster_data <- data.frame(
  run_name = rep(NA, data_length),
  cluster_ID = rep(NA, data_length),
  family_ID = rep(NA, data_length),
  seq_counts = rep(NA, data_length),
  read_counts = rep(NA, data_length),
  read_abun = rep(NA, data_length),
  sequence = rep(NA, data_length),
  cluster_color = rep(NA, data_length)
)

# initialize the index
index <- 1

# loop over each run
for (run_num in min(run_list):max(run_list)) {
  # loop over each cluster
  for (cluster_num in min(cluster_list):max(cluster_list)) {
    # set the current run, cluster, and family IDs
    cluster_data$run_name[index] <- run_num
    cluster_data$cluster_ID[index] <- cluster_num
    cluster_data$family_ID[index] <- fam_data[fam_data$cluster_ID == cluster_num, "family_ID"]
    # add sequence counts
    cluster_data$seq_counts[index] <-  nrow(seqs_identities[seqs_identities$peak_cluster_ID == cluster_num & seqs_identities$run_name == run_num,])
    # add read counts
    cluster_data$read_counts[index] <-  sum(seqs_identities[seqs_identities$peak_cluster_ID == cluster_num & seqs_identities$run_name == run_num,"read_counts"])
    # add read abundances
    cluster_data$read_abun[index] <- 100*cluster_data$read_counts[index]/quality[run_num]
    # add sequence
    cluster_data$sequence[index] <- fam_data[fam_data$cluster_ID == cluster_num, "sequence"]
    # add cluster color
    cluster_data$cluster_color[index] <- fam_data[fam_data$cluster_ID == cluster_num, "cluster_color"]
    # increment index
    index <- index+1
  }
}

# sort cluster data
cluster_data <- cluster_data[order(cluster_data$read_abun, decreasing = TRUE),]  

# line plot with percent abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_data, aes(x=as.character(run_name), y=read_abun, group=family_ID, color=cluster_color))+
  geom_line(size = 1) +
  geom_point() +
  theme_classic() +
  scale_color_identity(name = "Family", labels = cluster_data$family_ID, breaks = cluster_data$cluster_color, guide = "legend") +
  ylab("Percent Abundance") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/family_percent_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(cluster_abundances_plot)
dev.off()

# export plotting data
write.csv(cluster_data, file = paste(out_dir, "/family_abundances.csv", sep = ""), row.names = FALSE, quote = FALSE)
