#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N cluster_jobOutput

# script to cluster sequences using clustalo
# usage: qsub 08_cluster.sh 

# load the egapx software module (contains nextflow)
module load bio/2.0

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
clusterOut=$outputsPath"/clustered"

# move to the new directory
cd $clusterOut

# status message
echo "Beginning analysis..."

# filter to keep sequences with matching up- and down-stream sequences
clustalo -v -i $clusterOut"/aligned.fasta" -o $clusterOut"/clustered.fasta"

# status message
echo "Analysis complete!"
