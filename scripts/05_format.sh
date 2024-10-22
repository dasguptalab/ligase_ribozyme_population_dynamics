#!/bin/bash

# script to subset sequences and format headers
# usage: bash 05_format.sh analysisSubType
# usage: bash 05_format.sh a
# usage: bash 05_format.sh b
# usage: bash 05_format.sh c

# retrieve input analysis type
analysisSubType=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/cleaned_"$analysisSubType"_"$analysisTag

# make a new directory for analysis
filterOut=$outputsPath"/formatted_"$analysisSubType"_"$analysisTag
mkdir $formattedOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $formattedOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to outputs directory
cd $formattedOut

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
	cat $f1 | awk 'NR%2==0' | sort | uniq -c | sort -nrk1 > $formattedOut"/"$newName"_counts.tmp.txt"
	# for fastq files
	#cat $f1 | awk 'NR%4==2' | sort | uniq -c | sort -nrk1 > $formattedOut"/"$newName"_counts.tmp.txt"
	# get the run tag
	runTag=$(basename $f1 | cut -d"_" -f1)
	# get length of output file
	outLength=$(wc -l $formattedOut"/"$newName"_counts.tmp.txt" | cut -d"/" -f1 | tr -d " ")
	# make a file with the run tag on each line
	yes . | head -n $outLength | sed "s/./>$runTag/g" > $formattedOut"/"$newName"_run.tmp.txt"
	# make a file with a sequence of numbers to use as unique sequence IDs
	seq -f %1.0f 1 $outLength > $formattedOut"/"$newName"_ID.tmp.txt"
	# combine output counts and sequence file with run tags and convert to csv
	paste -d" " $formattedOut"/"$newName"_run.tmp.txt" $formattedOut"/"$newName"_ID.tmp.txt" $formattedOut"/"$newName"_counts.tmp.txt" | tr -s ' ' | sed "s/ /,/g" > $formattedOut"/"$newName".tmp.txt"
	# filter to remove sequences with less than 10 reads
	awk -F ',' '($3 > 10)' $formattedOut"/"$newName".tmp.txt" > $formattedOut"/"$newName"_above10.tmp.txt"
	# cut out the header data
	cut -d"," -f1-3 $formattedOut"/"$newName".tmp.txt" > $formattedOut"/"$newName"_header.tmp.txt"
	cut -d"," -f1-3 $formattedOut"/"$newName"_above10.tmp.txt" > $formattedOut"/"$newName"_above10_header.tmp.txt"
	# cut out the seqeunce data
	cut -d"," -f4 $formattedOut"/"$newName".tmp.txt" > $formattedOut"/"$newName"_seq.tmp.txt"
	cut -d"," -f4 $formattedOut"/"$newName"_above10.tmp.txt" > $formattedOut"/"$newName"_above10_seq.tmp.txt"
	# interleave the seqeunce data with the headers
	paste -d'\n' $formattedOut"/"$newName"_header.tmp.txt" $formattedOut"/"$newName"_seq.tmp.txt" > $formattedOut"/"$newName".fa"
	paste -d'\n' $formattedOut"/"$newName"_above10_header.tmp.txt" $formattedOut"/"$newName"_above10_seq.tmp.txt" > $formattedOut"/"$newName"_above10.fa"
	# clean up
	rm $formattedOut"/"$newName*".tmp."*
done

# status message
echo "Analysis complete!"
