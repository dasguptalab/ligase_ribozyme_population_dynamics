#!/bin/bash

# script to clean reads and keep only the variable 40bp region
# usage: bash 04_clean.sh analysisSubType
# usage: bash 04_clean.sh a
# usage: bash 04_clean.sh b
# usage: bash 04_clean.sh c

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# retrieve input analysis type
analysisSubType=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/filtered_"$analysisSubType"_"$analysisTag

# make a new directory for analysis
outputsPath=$outputsPath"/cleaned_"$analysisSubType"_"$analysisTag
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsPath

# status message
echo "Beginning analysis..."

# loop through all samples
for f1 in $inputsPath"/"*\.fq; do
	# status message
	echo "Processing $f1"
	# trim to sample tag
	newName=$(basename $f1 | sed 's/_filtered\.fq/_cleaned\.fa/')
	# filter to keep sequences with matching up- and down-stream sequences
	cat $f1 | sed "s/^.*GGACAGCG/START/g" | sed "s/CGCTGTCC.*$/END/g" | grep -Ex -B1 'START.{40}END' | grep -v "^--$" | sed "s/START//g" | sed "s/END//g" > $outputsPath"/"$newName
done

# status message
echo "Analysis complete!"
