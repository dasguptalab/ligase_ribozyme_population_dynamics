#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_jobOutput
#$ -pe smp 8

# script to cluster sequences using clustalo
# usage: qsub 07_cluster.sh sampleTag
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_above2\.fa; do echo $i; qsub 07_cluster.sh $i; done
## jobs 1541791 to 1541805

# load the software module
module load bio/0724

# set the sample file
sampleFile=$1

# retrieve the same tag
sampleTag=$(basename $sampleFile | sed 's/_formatted_above2\.fa//')
sampleTag=$sampleTag

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
outputsPath=$outputsPath"/07_clustered"
mkdir $outputsPath

# move to the new directory
cd $outputsPath

# status message
echo "Processing $sampleFile ..."

# cluster sequences
#clustalo --threads=$NSLOTS -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=500 
#clustalo -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=500 --full --percent-id --distmat-out=$outputsPath"/"$sampleTag"_distances.txt"
clustalo -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=1500

# status message
echo "Analysis complete!"
