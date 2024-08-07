#!/bin/bash

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: bash 06_analyzeReads.sh 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set inputs directory
formatOut=$outputsPath"/formatted"

# move to the new directory
cd $formatOut

# print read lengths
cat $formatOut"/combined.fmt.fasta" | awk 'NR % 2 == 0 {print length($0)}' | sort -n | uniq -c > $formatOut"/read_lengths.txt"

# print read counts
cat $formatOut"/combined.fmt.fasta" | awk 'NR % 2 == 0' | sort -n | uniq -c | sort -rk1 > $formatOut"/read_counts.txt"

# status message
echo "Analysis complete!"
