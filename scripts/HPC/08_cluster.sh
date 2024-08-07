#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_jobOutput
#$ -pe smp 32

# script to cluster sequences using clustalo
# usage: qsub 08_cluster.sh inputFile
# usage ex: qsub 08_cluster.sh combined.flt40.fmt.fasta
# usage ex: qsub 08_cluster.sh combined.flt.fmt.fasta

# load the software module
module load bio/0724

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set directory for inputs
formatOut=$outputsPath"/formatted"

# clean up input file name
nameTag=$(echo $inputFile | sed "s/\.fasta//g" | sed "s/\./_/g")

# make a new directory for analysis
clusterOut=$outputsPath"/clustered_"$nameTag
mkdir $clusterOut
# check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $clusterOut directory already exsists... please remove before proceeding."
#	exit 1
#fi

# move to the new directory
cd $clusterOut

# status message
echo "Beginning analysis of $nameTag ..."

# filter to keep sequences with matching up- and down-stream sequences
clustalo --threads=$NSLOTS -v -i $formatOut"/"$inputFile -o $clusterOut"/clustered_"$inputFile

# status message
echo "Analysis complete!"
