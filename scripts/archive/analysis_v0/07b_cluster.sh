#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_b_jobOutput
#$ -pe smp 8

# script to cluster sequences using clustalo
# usage: qsub 07b_cluster.sh

# load the software module
#module load bio/0724

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/06_formatted"

# make a new directory for analysis
outputsPath=$outputsPath"/07_clustered"
mkdir $outputsPath
outputsPath=$outputsPath"/07_clustered/07b_clustered"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsPath

# loop through all samples
for f1 in $inputsPath"/"*_above9\.fa; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_above9\.fa//')
	# status message
	echo "Processing $sampleTag ..."
	# cluster sequences
	#clustalo --threads=$NSLOTS -i $f1 --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" 
	#clustalo -i $f1 --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --full --percent-id --distmat-out=$outputsPath"/"$sampleTag"_distances.txt"
	clustalo -i $f1 --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa"
done

# status message
echo "Analysis complete!"
