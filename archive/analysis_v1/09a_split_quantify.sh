#!/bin/bash

# script to count the number of sequences shared across runs
# usage: bash 09a_split_quantify.sh inputRun runName
# usage ex: bash 09a_quantify.sh r7_S7_L001 r8_S8_L001

# retrieve input run name
inputRun=$1

# retrieve comparison run name
runName=$2

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve input sequences
inputSeqs=$outputsPath"/06_formatted/"$inputRun"_formatted.fa" ## quantification of all sequencess
#inputSeqs=$outputsPath"/06_formatted/"$inputRun"_formatted_above2.fa"

# name of a new directory for analysis
tablesOut=$outputsPath"/09a_quantified_all"

# make a new directory for analysis
mkdir $tablesOut

# setup tmp inputs data
inputSeqsData=$tablesOut"/"$inputRun"_formatted.tmp"
split --lines=10000 $inputSeqs $inputSeqsData".split.fa"

# initialize split counter
splitNum=0

# loop over each split file
for inputSeqs in $inputSeqsData".split.fa"*; do
	# increment split number
	splitNum=$(($splitNum+1))
	# submit job script
	qsub 09a_quantify.sh $inputRun $runName $inputSeqs $splitNum
done
