#!/bin/bash

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: bash 06_analyzeReads.sh 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# set trimming directory
clusterOut=$outputsPath"/clustered"

# move to the new directory
cd $clusterOut

# print read lengths
cat $clusterOut"/combined.fasta" | awk 'NR % 2 == 0 {print length($0)}' | sort -n | uniq -c > $clusterOut"/read_lengths.txt"

# print read counts
cat $clusterOut"/combined.fasta" | awk 'NR % 2 == 0' | sort -n | uniq -c | sort -rk1 > $clusterOut"/read_counts.txt"

# status message
echo "Analysis complete!"
