#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster40_jobOutput
#$ -pe smp 32

# script to cluster sequences using clustalo
# usage: qsub 08a_cluster40.sh 

# load the software module
#module load bio
module load bio/0724

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set directory for inputs
formatOut=$outputsPath"/formatted"

# make a new directory for analysis
clusterOut=$outputsPath"/clustered40"
mkdir $clusterOut
# check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $clusterOut directory already exsists... please remove before proceeding."
#	exit 1
#fi

# move to the new directory
cd $clusterOut

# status message
echo "Beginning analysis of $formatOut/combined.flt40.fmt.fasta ..."

# filter to keep sequences with matching up- and down-stream sequences
#clustalo --threads=32 -v -i $formatOut"/combined.flt40.fmt.fasta" -o $clusterOut"/clustered.flt40.fmt.fasta"
clustalo --threads=$NSLOTS -v -i $formatOut"/combined.flt40.fmt.fasta" -o $clusterOut"/clustered.flt40.fmt.fasta"

# status message
echo "Analysis complete!"
