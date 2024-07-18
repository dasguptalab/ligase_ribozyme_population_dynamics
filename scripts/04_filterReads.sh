#!/bin/bash

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: bash 04_filterReads.sh 

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# set trimming directory
trimOut=$outputsPath"/trimmed"

# make a new directory for analysis
filterOut=$outputsPath"/filtered"
mkdir $filterOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $filterOut directory already exsists... please remove before proceeding."
	exit 1
fi

# unzip all paired reads
gunzip $trimOut"/"*_p* $trimOut

# move to the new directory
cd $filterOut

# status message
echo "Beginning analysis..."

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in $trimOut"/"*_p*\.fq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/\.fq//')
	# status message
	echo "Processing $sampleTag"
	# filter to keep sequences with matching up- and down-stream sequences
	cat $f1 | grep -E -B1 -A2 "^GGACAGCG.*CGCTGTCC.*" | grep -v "^--$" | sed "s/^GGACAGCG//g" | sed "s/CGCTGTCC.*$//g" > $filterOut"/"$sampleTag".flt.fq"
	# status message
	echo "$sampleTag processed!"
done

# status message
echo "Analysis complete!"
