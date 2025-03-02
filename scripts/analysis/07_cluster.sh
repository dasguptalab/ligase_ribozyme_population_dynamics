#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_jobOutput
#$ -pe smp 8

# script to cluster sequences using clustalo
# usage: qsub 07_cluster.sh sampleTag
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_formatted\.fa; do sampleTag=$(basename $i | sed 's/_formatted\.fa//'); echo $sampleTag; qsub 07_cluster.sh $sampleTag; done
## jobs 

# load the software module
module load bio/0724

# retrieve the same tag
sampleTag=$1

# set the sample file
sampleFile=$1"_formatted.fa"

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/06_formatted"

# make a new directory for analysis
outputsPath=$outputsPath"/07_clustered"
mkdir $outputsPath

# move to the new directory
cd $outputsPath

# status message
echo "Processing $sampleTag ..."

# cluster sequences
#clustalo --threads=$NSLOTS -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=500 
#clustalo -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=500 --full --percent-id --distmat-out=$outputsPath"/"$sampleTag"_distances.txt"
clustalo -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=500

# status message
echo "Analysis complete!"
