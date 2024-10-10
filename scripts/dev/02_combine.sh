#!/bin/bash

# script to combine files of merged trimmed paired reads with trimmed unpaired reads
# usage: bash 02_combine.sh analysisType
# usage: bash 02_combine.sh trimmed_s4q20

# retrieve input analysis type
analysisType=$1

# retrieve analysis tag
analysisTag=$(echo $1 | sed "s/trimmed_//g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set the inputs directory
inputsPath=$outputsPath

# make a new directory for analysis
outputsCombined=$outputsPath"/combined_"$analysisTag
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
gunzip -v $inputsPath"/trimmed_"$analysisTag"/"*\.gz
gunzip -v $inputsPath"/merged_"$analysisTag"/"*\.gz

# loop over un-filtered merged reads for each run
for f1 in $inputsPath"/trimmed_"$analysisTag"/"*_u*\.fq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_stiched_reads\.fq//')
	# status message
	echo "Processing $sampleTag ..."
	# combine un-filtered merged,.fqiled merged, and unpaired trimmed reads
	cat $f1 $inputsPath"/trimmed_"$analysisTag"/"$sampleTag"_u"*\.fq $inputsPath"/merged_"$analysisTag"/"$sampleTag*\.fq >> $outputsCombined"/"$sampleTag"_combined.fq"
done

# double check that there are no duplicate reads
#for i in $outputsCombined"/"*_combined\.fq; do echo $i; cat $i | awk 'NR%2==1' | cut -d' ' -f1 | sort | uniq -c | sort -n | head; done

# status message
echo "Analysis complete!"
