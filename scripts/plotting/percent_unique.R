#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/07a_clustered"

# round numbers
rounds <- c(1, 2, 3, 4, 5, 6, 7, 8)

# % diversity per round
#diversity_doped <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20, 97.30, 92.40, 86.43)
diversity <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20)

# calculate % unique per round
unique_reads <- diversity

# create data frame for plotting
unique_rounds <- data.frame(
  round_nums = rounds,
  #round_nums = as.character(rounds),
  unique_percent = unique_reads
)

# line plot with % of unique sequences per round
unique_rounds_plot <- ggplot(data=unique_rounds, aes(round_nums, unique_percent)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylab("Percent Unique Sequences") +
  xlab("Round Number")
# save the plot
exportFile <- paste(out_dir, "/percent_unique.png", sep = "")
png(exportFile, units="in", width=5, height=5, res=300)
print(unique_rounds_plot)
dev.off()

# export plotting data
write.csv(unique_rounds, file = paste(out_dir, "/data/percent_unique.csv", sep = ""), row.names = FALSE, quote = FALSE)
