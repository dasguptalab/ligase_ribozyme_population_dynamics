#!/bin/bash

# script to combine files of merged trimmed paired reads with trimmed unpaired reads
# usage: bash 05_combine.sh 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set the inputs directory
inputsPath=$outputsPath

# status message
echo "Beginning analysis..."

# make a new directory for analysis
outputsFiltered=$outputsPath"/filtered_combined"
mkdir $outputsFiltered
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsFiltered directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsFiltered

# status message
echo "Analyzing filtered combined data..."

# loop over filtered merged reads for each run
for f1 in $inputsPath"/filtered/"*_uForward\.flt\.fq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_uForward\.flt\.fq//')
	# status message
	echo "Processing $sampleTag"
	# combine filtered merged,.fqiled merged, and unpaired trimmed reads
	cat $inputsPath"/filtered/"$sampleTag*\.fq >> $outputsFiltered"/"$sampleTag"_combined.fq"
	# status message
	echo "$sampleTag processed!"
done

# double check that there are no duplicate reads
#for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/combined_filtered/*_combined\.fq; do echo $i; cat $i | awk 'NR%2==1' | cut -d' ' -f1 | sort | uniq -c | sort -n | head; done

# make a new directory for analysis
outputsCombined=$outputsPath"/combined"
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

# loop over un-filtered merged reads for each run
for f1 in $inputsPath"/trimmed/"*_uForward\.fq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_uForward\.fq//')
	# status message
	echo "Processing $sampleTag"
	# combine un-filtered merged,.fqiled merged, and unpaired trimmed reads
	cat $f1 $inputsPath"/merged/"$sampleTag\_stiched_reads\.fq $inputsPath"/merged/logs/"$sampleTag\_stiched_reads.fqiled\.fq_*\.fastq >> $outputsCombined"/"$sampleTag"_combined.fq"
	# status message
	echo "$sampleTag processed!"
done

# double check that there are no duplicate reads
#for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/combined/*_combined\.fq; do echo $i; cat $i | awk 'NR%2==1' | cut -d' ' -f1 | sort | uniq -c | sort -n | head; done

# status message
echo "Analysis complete!"
