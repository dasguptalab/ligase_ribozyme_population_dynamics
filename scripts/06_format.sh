#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_subset_format_jobOutput
#$ -q largemem

# script to subset sequences and format headers
# usage: bash 06_format.sh inputFile

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set trimming directory
inputsPath=$outputsPath"/filtered_combined"

# make a new directory for analysis
formatOut=$outputsPath"/formatted_subset"
mkdir $formatOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $formatOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $formatOut

# create sequence run name
#newName=$(basename $f1 | sed "s/_L001_p.*\.flt\.fq//g")
# clean up file name
fileName=$(basename $inputFile | sed "s/\.fa/\.fmt\.fa/g")

# pre-clean up
rm $formatOut"/"$fileName

# status message
echo "Beginning analysis..."

# subset read sequences and re-format headers
while read -r line; do
	# check if the current line is a header
	if [[ $line == "@"* ]]; then # read header
		# save previous line, which should be the leader of the current sequence
		headerLine=$line
	else # read sequence
		# status message
		echo "Processing $line ..."
		# get read count
		readCount=$(cat $inputFile | grep $line | wc -l | tr -d " ")
		# check if there are more than 10 of the current read
		if (( $readCount > 10 )); then # add read info to outputs
			# re-format read headers and add to outputs
			echo $headerLine | cut -d" " -f1 | sed "s/^@/>/g" >> $formatOut"/"$fileName
			# add sequence to outputs
			echo $line >> $formatOut"/"$fileName
		fi
	fi
done < $inputFile

# status message
echo "Analysis complete!"
