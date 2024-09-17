#!/bin/bash

# script to combine files of merged trimmed paired reads with trimmed unpaired reads
# usage: bash 04_combine.sh analysisType

# retrieve input analysis type
analysisType=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set the inputs directory
inputsPath=$outputsPath

# retrieve the inputs path
inputsPath=$outputsPath"/"$analysisType

# make a new directory for analysis
outputsCombined=$outputsPath"/combined_"$analysisType
mkdir $outputsCombined
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsCombined directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsCombined

# status message
echo "Analyzing combined data..."

# unzip any gz read files
gunzip -v $inputsPath"/trimmed/"*\.gz
gunzip -v $inputsPath"/merged/"*\.gz
gunzip -v $inputsPath"/merged/logs/"*\.gz

# loop over un-filtered merged reads for each run
for f1 in $inputsPath"/"*_uForward\.fq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_uForward\.fq//')
	# status message
	echo "Processing $sampleTag ..."
	# combine un-filtered merged,.fqiled merged, and unpaired trimmed reads
	cat $f1 $inputsPath"/trimmed/"$sampleTag"_uReverse.fq" $inputsPath"/merged/"$sampleTag\_stiched_reads\.fq $inputsPath"/merged/logs/"$sampleTag\_stiched_reads_failed\.fq_*\.fastq >> $outputsCombined"/"$sampleTag"_combined.fq"
done

# double check that there are no duplicate reads
for i in $outputsPath"/combined/"*_combined\.fq; do echo $i; cat $i | awk 'NR%2==1' | cut -d' ' -f1 | sort | uniq -c | sort -n | head; done

# status message
echo "Analysis complete!"
