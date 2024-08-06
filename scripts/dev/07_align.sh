#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N align_jobOutput
#$ -pe smp 32

# script to align sequences using muscle
# usage: qsub 07_align.sh

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
alignOut=$outputsPath"/clustered"

# move to the new directory
cd $alignOut

# status message
echo "Beginning analysis..."

# align reads
muscle -threads 32 -super5 $alignOut"/combined.fasta" -output $alignOut"/aligned.fasta"
#muscle -in $alignOut"/combined.fasta" -out $alignOut"/aligned.fasta"

# status message
echo "Analysis complete!"
