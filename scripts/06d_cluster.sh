#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_maxnumseq_clustalo_jobOutput
#$ -pe smp 63
#$ -q largemem

# script to cluster sequences using clustalo and --maxnumseq=1300000
# usage: qsub 06d_cluster.sh inputFile
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06d_cluster.sh "${fileList[$i]}"; done

# retrieve input file
inputFile=$1

# retrieve software path
softwarePath=$(grep "software_cdhit:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/software_cdhit://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# clean up input file name
nameTag=$(basename $inputFile | sed "s/\.fa//g" | sed "s/\./_/g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# make a directory for the clustering outputs
clusterOut=$outputsPath"/clustered_cd_hit_"$analysisTag
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
$softwarePath"/"cd-hit -T $NSLOTS -sc 1 -sf 1 -i $inputFile -o $clusterOut"/clustered_"$nameTag

# status message
echo "Analysis complete!"
