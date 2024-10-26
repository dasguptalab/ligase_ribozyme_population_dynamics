#!/bin/bash

# script to subset sequences and format headers
# usage: bash 06_format.sh

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/combined"

# make a new directory for analysis
outputsPath=$outputsPath"/formatted"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to outputs directory
cd $outputsPath

# status message
echo "Beginning analysis..."

# loop through all samples
for f1 in $inputsPath"/"*\.fa; do
	# status message
	echo "Processing file: $f1"
	# trim to sample tag
	newName=$(basename $f1 | sed 's/_cleaned\.fa/_formatted/')
	# print read counts
	# for fasta files
	cat $f1 | awk 'NR%2==0' | sort | uniq -c | sort -nrk1 > $outputsPath"/"$newName"_counts.tmp.txt"
	# for fastq files
	#cat $f1 | awk 'NR%4==2' | sort | uniq -c | sort -nrk1 > $outputsPath"/"$newName"_counts.tmp.txt"
	# get the run tag
	runTag=$(basename $f1 | cut -d"_" -f1)
	# get length of output file
	outLength=$(wc -l $outputsPath"/"$newName"_counts.tmp.txt" | cut -d"/" -f1 | tr -d " ")
	# make a file with the run tag on each line
	yes . | head -n $outLength | sed "s/./>$runTag/g" > $outputsPath"/"$newName"_run.tmp.txt"
	# make a file with a sequence of numbers to use as unique sequence IDs
	seq -f %1.0f 1 $outLength > $outputsPath"/"$newName"_ID.tmp.txt"
	# combine output counts and sequence file with run tags and convert to csv
	paste -d" " $outputsPath"/"$newName"_run.tmp.txt" $outputsPath"/"$newName"_ID.tmp.txt" $outputsPath"/"$newName"_counts.tmp.txt" | tr -s ' ' | sed "s/ /,/g" > $outputsPath"/"$newName".tmp.txt"
	# filter to remove sequences with less than 10 reads
	awk -F ',' '($3 > 10)' $outputsPath"/"$newName".tmp.txt" > $outputsPath"/"$newName"_above10.tmp.txt"
	# cut out the header data
	cut -d"," -f1-3 $outputsPath"/"$newName".tmp.txt" > $outputsPath"/"$newName"_header.tmp.txt"
	cut -d"," -f1-3 $outputsPath"/"$newName"_above10.tmp.txt" > $outputsPath"/"$newName"_above10_header.tmp.txt"
	# cut out the seqeunce data
	cut -d"," -f4 $outputsPath"/"$newName".tmp.txt" > $outputsPath"/"$newName"_seq.tmp.txt"
	cut -d"," -f4 $outputsPath"/"$newName"_above10.tmp.txt" > $outputsPath"/"$newName"_above10_seq.tmp.txt"
	# interleave the seqeunce data with the headers
	paste -d'\n' $outputsPath"/"$newName"_header.tmp.txt" $outputsPath"/"$newName"_seq.tmp.txt" > $outputsPath"/"$newName".fa"
	paste -d'\n' $outputsPath"/"$newName"_above10_header.tmp.txt" $outputsPath"/"$newName"_above10_seq.tmp.txt" > $outputsPath"/"$newName"_above10.fa"
	# clean up
	rm $outputsPath"/"$newName*".tmp."*
done

# status message
echo "Analysis complete!"
