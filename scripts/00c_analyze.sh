#!/bin/bash

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: bash 00c_analyze.sh analysisType
# usage: bash 00c_analyze.sh trimmed
# usage: bash 00c_analyze.sh merged
# usage: bash 00c_analyze.sh filtered
# usage: bash 00c_analyze.sh combined
# usage: bash 00c_analyze.sh combined_filtered

# retrieve input analysis type
analysisType=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# check input analysis type
if [[ $analysisType == "raw" ]]; then
	# retrieve raw paired reads absolute path for alignment
	readPath=$(grep "pairedReads:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/pairedReads://g")
else 
	# set reads path
	readPath=$outputsPath"/"$analysisType
fi

#Make a new directory for analysis
analysisOut=$outputsPath"/analyzed_"$analysisType
mkdir $analysisOut
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $analysisOut directory already exsists... please remove before proceeding."
	exit 1
fi

#Move to the new directory
cd $analysisOut

# status message
echo "Beginning analysis..."

# name output files
lengthsOut=$analysisOut"/read_lengths.txt"
countsOut=$analysisOut"/read_counts.txt"

# add headers
echo "run,count,length" > $lengthsOut
echo "run,count,sequence" > $countsOut
#echo "run,total,quality,unique,diversity"

# loop over each file
for f1 in $readPath"/"*; do
	# status message
	echo "Processing file: $f1"
	# check analysis type
	if [[ $analysisType == *"filtered" ]]; then
		# print read lengths
		cat $f1 | awk 'NR%2==0 {print length($0)}' | sort -n | uniq -c > $lengthsOut".tmp.txt"
		# print read counts
		cat $f1 | awk 'NR%2==0' | sort -n | uniq -c | sort -rk1 > $countsOut".tmp.txt"
	else
		# print read lengths
		cat $f1 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c > $lengthsOut".tmp.txt"
		# print read counts
		cat $f1 | awk 'NR%4==2' | sort -n | uniq -c | sort -rk1 > $countsOut".tmp.txt"
	fi
	# get the run tag
	runTag=$(basename $f1 | cut -d"_" -f1)
	# get length of output file
	outLength=$(wc -l $lengthsOut".tmp.txt" | cut -d"/" -f1 | tr -d " ")
	# make a file with the run tag on each line
	yes . | head -n $outLength | sed "s/./$runTag/g" > $lengthsOut".run.tmp.txt"
	# combine output file with run tags and convert to csv
	paste -d" " $lengthsOut".run.tmp.txt" $lengthsOut".tmp.txt" | tr -s ' ' | sed "s/ /,/g" >> $lengthsOut
	# get length of output file
	outLength=$(wc -l $countsOut".tmp.txt" | cut -d"/" -f1 | tr -d " ")
	# make a file with the run tag on each line
	yes . | head -n $outLength | sed "s/./$runTag/g" > $countsOut".run.tmp.txt"
	# combine output file with run tags and convert to csv
	paste -d" " $countsOut".run.tmp.txt" $countsOut".tmp.txt" | tr -s ' ' | sed "s/ /,/g" >> $countsOut
	# clean up
	rm $lengthsOut".tmp.txt"
	rm $lengthsOut".run.tmp.txt"
	rm $countsOut".tmp.txt"
	rm $countsOut".run.tmp.txt"
	# status message
	echo "File processed!"
done

# status message
echo "Analysis complete!"
