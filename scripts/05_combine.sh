#!/bin/bash

# script to combine files of merged trimmed paired reads with trimmed unpaired reads
# usage: bash 02_combine.sh

# retrieve analysis tag
analysisTag=$(echo $1 | sed "s/trimmed_//g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/cleaned"

# make a new directory for analysis
outputsPath=$outputsPath"/combined"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsCombined

# unzip any gz read files
gunzip -v $inputsPath"/"*

# loop through all samples
for f1 in $inputsPath"/"*_trimmed\.fa; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_trimmed\.fa//')
	# status message
	echo "Processing $sampleTag ..."
	# combine un-filtered merged,.fqiled merged, and unpaired trimmed reads
	cat $f1 $inputsPath"/"$sampleTag"_pForward.fq" $inputsPath"/"$sampleTag"_uForward.fq" >> $outputsPath"/"$sampleTag"_combined.fa"
done

# double check that there are no duplicate reads
#for i in $outputsCombined"/"*_combined\.fa; do echo $i; cat $i | awk 'NR%2==1' | cut -d' ' -f1 | sort | uniq -c | sort -n | head; done

# status message
echo "Analysis complete!"
