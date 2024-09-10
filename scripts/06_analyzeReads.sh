#!/bin/bash

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: bash 06_analyzeReads.sh 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set inputs directory
filterOut=$outputsPath"/filteredByPair"

# move to the new directory
cd $filterOut

# name output files
lengthsOut=$filterOut"/read_lengths.txt"
countsOut=$filterOut"/read_counts.txt"

# add headers
echo "run,count,length" > $lengthsOut
echo "run,count,sequence" > $countsOut
echo "run,total,quality,unique,diversity"

# loop over each file
for f1 in $filterOut"/"*; do
	# print read lengths
	cat $f1 | awk 'NR % 2 == 0 {print length($0)}' | sort -n | uniq -c > $lengthsOut".tmp.txt"
	# get length of output file
	outLength=$(wc -l $lengthsOut".tmp.txt")
	# make a file with the run tag on each line
	yes . | head -n $outLength > $lengthsOut".run.tmp.txt"
	# combine output file with run tags and convert to csv
	paste -d" " $lengthsOut".run.tmp.txt" $lengthsOut".tmp.txt" | sed "s/ /,/g" >> $lengthsOut
	# print read counts
	cat $f1 | awk 'NR % 2 == 0' | sort -n | uniq -c | sort -rk1 > $countsOut".tmp.txt"
	# get length of output file
	outLength=$(wc -l $countsOut".tmp.txt")
	# make a file with the run tag on each line
	yes . | head -n $outLength > $countsOut".run.tmp.txt"
	# combine output file with run tags and convert to csv
	paste -d" " $countsOut".run.tmp.txt" $countsOut".tmp.txt" | sed "s/ /,/g" >> $countsOut
	# clean up
	rm $lengthsOut".tmp.txt"
	rm $lengthsOut".run.tmp.txt"
	rm $countsOut".tmp.txt"
	rm $countsOut".run.tmp.txt"
done

# status message
echo "Analysis complete!"
