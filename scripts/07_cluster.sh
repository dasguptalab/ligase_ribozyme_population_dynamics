#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_jobOutput
#$ -pe smp 8
#$ -q largemem

# script to cluster sequences using clustalo
# usage: qsub 07_cluster.sh inputFile
# usage ex: qsub 07_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs_flash/formatted/r8_S8_L001_formatted.fa
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs_flash/formatted/*_above9); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 07_cluster.sh "${fileList[$i]}"; done
# alternate usage: bash 07_cluster.sh inputFile
# usage ex: bash 07_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs_flash/formatted/r8_S8_L001_formatted_above9.fa
# usage ex: bash 07_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs_flash/formatted/r7_S7_L001_formatted_above9.fa
# usage ex: bash 07_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs_flash/formatted/r6_S6_L001_formatted_above9.fa

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
clusterOut=$outputsPath"/clustered"
#clusterOut=$outputsPath"/clustered_trimmed"
mkdir $clusterOut

# make a new directory for analysis
clusterOut=$clusterOut"/"$nameTag
mkdir $clusterOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $clusterOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $clusterOut

# status message
echo "Beginning analysis of $nameTag ..."

# cluster sequences
#clustalo --threads=$NSLOTS -i $inputFile --clustering-out=$clusterOut"/"$nameTag"_clustered.aux" -o $clusterOut"/"$nameTag"_aligned.fa" --cluster-size=500
clustalo --threads=8 -i $inputFile --clustering-out=$clusterOut"/"$nameTag"_clustered.aux" -o $clusterOut"/"$nameTag"_aligned.fa" --cluster-size=500

# status message
echo "Analysis complete!"
