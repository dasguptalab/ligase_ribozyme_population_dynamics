#!/bin/bash

# script to create files with randomized nucelotide sequences
# usage: bash 14_randomize.sh

# retrieve analysis outputs absolute path
outputsPath="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs"
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/13_randomized_sequences"

# create outputs directory
mkdir $outDir

# set of nucleotides 
nucSet="ACTG"

# length of seqeuences
seqLength=40

# length of file
numSeqs=2000000

# name out output file
outFile=$outDir"/randomized_sequences.fa"

# add the header
echo > $outFile

# loop over each line of the file
for seqNum in $(seq 1 $numSeqs); do
	# status message
	echo "Processing sequence number $seqNum"
	# reset the current sequence
	currSeq=""
	# loop over each nucleotide
	for nucNum in $(seq 1 $seqLength); do
		# generate a random number
		nucIndex=$((RANDOM % 4))
		# select a random nucleotide
		randomNuc="${nucSet:$nucIndex:1}"
		# add the random nucelotide to the sequence
		currSeq="$currSeq$randomNuc"
	done
	# add the randomized sequence to the file
	echo $currSeq >> $outFile
done
