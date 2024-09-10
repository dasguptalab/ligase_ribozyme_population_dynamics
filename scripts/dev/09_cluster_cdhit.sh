#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_cdhit_jobOutput

# script to cluster sequences using cd-hit
# usage: qsub 09_cluster_cdhit.sh inputFile
# usage ex: qsub 09_cluster_cdhit.sh combined_noDoped.flt.fmt.fasta
## job 
# usage ex: qsub 09_cluster_cdhit.sh combined_noDoped.flt40.fmt.fasta
## job 750090 -> run1

# retrieve input file
inputFile=$1

# retrieve software path
softwarePath=$(grep "software_cdhit:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/software_cdhit://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set directory for inputs
formatOut=$outputsPath"/formatted"

# clean up input file name
nameTag=$(echo $inputFile | sed "s/\.fasta//g" | sed "s/\./_/g")

# make a new directory for analysis
clusterOut=$outputsPath"/clustered_cdhit_"$nameTag
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

# cluster sequences by similarity
$softwarePath"/"cd-hit-para.pl --Q 64 --T "SGE" -i $formatOut"/"$inputFile -o $clusterOut"/clustered_"$inputFile --R $clusterOut"/restart.log"

# status message
echo "Analysis complete!"
