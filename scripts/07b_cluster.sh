#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_b_jobOutput
#$ -pe smp 8
#$ -q largemem

# script to cluster sequences using clustalo
# usage: qsub 07b_cluster.sh inputFile
# usage ex: qsub 07b_cluster.sh r8_S8_L001_formatted.fa
# alternate usage: bash 07b_cluster.sh inputFile
# usage ex: bash 07b_cluster.sh r8_S8_L001_formatted_above9.fa

# load the software module
module load bio/0724

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# clean up input file name
nameTag=$(basename $inputFile | sed "s/\.fa//g" | sed "s/\./_/g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# make a directory for the clustering outputs
clusterOut=$outputsPath"/07_clustered"
mkdir $clusterOut

# move to the new directory
cd $clusterOut

# status message
echo "Beginning analysis of $nameTag ..."

# cluster sequences
#clustalo --threads=$NSLOTS -i $inputFile --clustering-out=$clusterOut"/"$nameTag"_clustered.aux" -o $clusterOut"/"$nameTag"_aligned.fa"
clustalo --threads=8 -i $inputFile --clustering-out=$clusterOut"/"$nameTag"_clustered.aux" -o $clusterOut"/"$nameTag"_aligned.fa"

# status message
echo "Analysis complete!"
