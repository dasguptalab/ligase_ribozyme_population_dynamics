
# load libraries
library(ggplot2)

# import length data
lengthData <- read.delim("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/clustered/read_lengths.fmt.txt", sep=" ")

# plot lengths scatter
ggplot(data=lengthData, mapping = aes(y=count, x=length)) +
  geom_line() +
  theme_minimal()

# plot lengths histogram
ggplot(data=lengthData, mapping = aes(x=length)) +
  geom_histogram(aes(y = after_stat(count))) +
  theme_minimal()
