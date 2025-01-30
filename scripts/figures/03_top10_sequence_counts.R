#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
library(dplyr)
library(ComplexHeatmap)

# supress cowplot package messages
suppressMessages( require(cowplot) )

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/03_top10_sequence_counts"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# round numbers
rounds <- c(1, 2, 3, 4, 5, 6, 7, 8)

# % diversity per round
#diversity_doped <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20, 97.30, 92.40, 86.43)
#diversity <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20)

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
#quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in sequence count data
seqs_counts <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11_quantified_top10/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))

# reverse complement the sequences
#seqs_counts$sequence <- rev(chartr("ATGC","TACG",seqs_counts$sequence))

# add the fraction abundances
#seqs_counts$frac_abundance <- NA
#for (run_num in 1:8) {
#  seqs_counts[seqs_counts$counts_run_name == run_num, "frac_abundance"] <- seqs_counts[seqs_counts$counts_run_name == run_num, "counts"]/quality[run_num]
#}

# add log values
seqs_counts$log_counts <- log(seqs_counts$counts)
#seqs_counts$log_frac_abundance <- log(seqs_counts$frac_abundance)

# add percent abundance
#seqs_counts$perc_abundance <- 100*seqs_counts$frac_abundance

# set infinite and NA values equal to zero
is.na(seqs_counts)<-sapply(seqs_counts, is.infinite)

# re-order the data for plotting
#seqs_counts <- seqs_counts[order(seqs_counts$log_counts, decreasing=TRUE),]

# setup midpoint values for plotting
seqs_counts_noNA <- seqs_counts
seqs_counts_noNA[is.na(seqs_counts_noNA)] <- 0
mid_log_counts <- max(seqs_counts_noNA$log_counts)/2
#mid_log_frac_abundance <- log(min(seqs_counts_noNA[seqs_counts_noNA$frac_abundance != 0,"frac_abundance"]))/2
#mid_perc_abundance <- max(seqs_counts$perc_abundance)/2

# create ID factor list
ID_list <- rep("1",8)
#ID_list <- rep(1,8)
for (i in 2:10) {
  ID_list <- c(ID_list, rep(i,8))
}
ID_data <- as.factor(ID_list)
#ID_data <- factor(ID_data, levels = levels(ID_data)[c(1, 3:10, 2)])
ID_data <- factor(ID_data, levels = levels(ID_data)[rev(c(1, 3:10, 2))])

# initialize plot data
counts_heatmap_list <- vector('list', 8)

# loop over each round and create heatmaps
for (run_num in 1:8) {
  # subset seq data
  seqs_counts_subset <- seqs_counts[seqs_counts$run_name == run_num,]
  # remove run tags
  seqs_counts_subset[c("run", "ID")] <- do.call(rbind, strsplit(seqs_counts_subset$sequence_ID, "_"))
  #seqs_counts_subset[c("run", "Sequence")] <- do.call(rbind, strsplit(seqs_counts_subset$sequence, "_"))
  # add sequence top 10 ID
  seqs_counts_subset$ID <- ID_data
  # set round plot title
  run_title <- paste("Round", run_num)
  # To-do: consider using factors
  # create heatmap of log counts
  #counts_heatmap_subset <- ggplot(data = seqs_counts_subset, aes(counts_run_name, reorder(as.character(ID), log_counts), fill= log_counts)) + 
  counts_heatmap_subset <- ggplot(data = seqs_counts_subset, aes(counts_run_name, ID, fill= log_counts)) + 
  #counts_heatmap_subset <- ggplot(data = seqs_counts_subset, aes(counts_run_name, reorder(as.character(Sequence), log_counts), fill= log_counts)) + 
    theme_bw() +
    geom_tile(colour = "black") +
    ggtitle(run_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    #ylab("Sequence ID") +
    ylab("Sequence") +
    xlab("Round Number") +
    scale_fill_gradient2(name = "Log Counts",
                         low = safe_colors[3],
                         mid = safe_colors[4],
                         high = safe_colors[5],
                         midpoint = mid_log_counts,
                         na.value = "white") #+
    #scale_y_discrete(labels=seqs_counts_subset$ID)
  # save the plot
  exportFile <- paste(out_dir, "/r", run_num, "_top10_sequence_log_counts.png", sep = "")
  #png(exportFile, units="in", width=5, height=5, res=300)
  png(exportFile, units="in", width=10, height=5, res=300)
  print(counts_heatmap_subset)
  dev.off()
  # add the log counts heatmap to the plot list
  counts_heatmap_list[[run_num]] <- counts_heatmap_subset
}

# combine heatmaps with the log counts for each of the top 10 sequences per round
grid_log_counts_plot <- plot_grid(plotlist=counts_heatmap_list,  ncol = 4, align = 'v')
# save the plot
exportFile <- paste(out_dir, "/top10_sequence_log_counts.png", sep = "")
png(exportFile, units="in", width=40, height=10, res=300)
print(grid_log_counts_plot)
dev.off()

# export plotting data
write.csv(seqs_counts, file = paste(out_dir, "/top10_sequence_counts.csv", sep = ""), row.names = FALSE, quote = FALSE)

# create table of top 10 sequences per run
sequence_data <- data.frame(
  run_1 = seqs_counts[seqs_counts$run_name == 1 & seqs_counts$counts_run_name == 8, "sequence"],
  run_2 = seqs_counts[seqs_counts$run_name == 2 & seqs_counts$counts_run_name == 8, "sequence"],
  run_3 = seqs_counts[seqs_counts$run_name == 3 & seqs_counts$counts_run_name == 8, "sequence"],
  run_4 = seqs_counts[seqs_counts$run_name == 4 & seqs_counts$counts_run_name == 8, "sequence"],
  run_5 = seqs_counts[seqs_counts$run_name == 5 & seqs_counts$counts_run_name == 8, "sequence"],
  run_6 = seqs_counts[seqs_counts$run_name == 6 & seqs_counts$counts_run_name == 8, "sequence"],
  run_7 = seqs_counts[seqs_counts$run_name == 7 & seqs_counts$counts_run_name == 8, "sequence"],
  run_8 = seqs_counts[seqs_counts$run_name == 8 & seqs_counts$counts_run_name == 8, "sequence"]
)

# export sequence data
write.csv(sequence_data, file = paste(out_dir, "/top10_sequences.csv", sep = ""), row.names = FALSE, quote = FALSE)

# setup and create upset plot
lt = list(
          run_1 = seqs_counts[seqs_counts$run_name == 1 & seqs_counts$counts_run_name == 8, "sequence"],
          run_2 = seqs_counts[seqs_counts$run_name == 2 & seqs_counts$counts_run_name == 8, "sequence"],
          run_3 = seqs_counts[seqs_counts$run_name == 3 & seqs_counts$counts_run_name == 8, "sequence"],
          run_4 = seqs_counts[seqs_counts$run_name == 4 & seqs_counts$counts_run_name == 8, "sequence"],
          run_5 = seqs_counts[seqs_counts$run_name == 5 & seqs_counts$counts_run_name == 8, "sequence"],
          run_6 = seqs_counts[seqs_counts$run_name == 6 & seqs_counts$counts_run_name == 8, "sequence"],
          run_7 = seqs_counts[seqs_counts$run_name == 7 & seqs_counts$counts_run_name == 8, "sequence"],
          run_8 = seqs_counts[seqs_counts$run_name == 8 & seqs_counts$counts_run_name == 8, "sequence"]
)
m = make_comb_mat(lt)
cs = comb_size(m)
ht = UpSet(m, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF", top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))))
ht = draw(ht)
# save the plot
exportFile <- paste(out_dir, "/top10_sequences_upset.png", sep = "")
png(exportFile, units="in", width=40, height=10, res=300)
print(ht)
dev.off()
