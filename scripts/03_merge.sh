#!/bin/bash

# script to perform merging of paired end reads into single reads
# usage: bash 03_merge.sh

# retrieve input analysis type
analysisType=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

#Make a new directory for analysis
qcOut=$outputsPath"/qc_"$analysisType
mkdir $qcOut
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $qcOut directory already exsists... please remove before proceeding."
	exit 1
fi

#Move to the new directory
cd $qcOut

# status message
echo "Processing..."

# perform merging
#fastqc $readPath"/"* -o $qcOut
NGmerge -v -n 8 -1 <file> -2 <file> -o <file> -o stiched_reads.fa -m 20 -p 0 -l log_stitching_results.txt -f stiched_reads_failed.fa -j log_formatted_alignments.txt -q 33 -u 40 -t ' '

#Print status message
echo "Analysis complete!"
