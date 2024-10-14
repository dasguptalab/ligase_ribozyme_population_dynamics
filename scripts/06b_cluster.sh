#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_500_clustalo_jobOutput
#$ -pe smp 63
#$ -q largemem

# script to cluster sequences using clustalo and --cluster-size=500
# usage: qsub 06b_cluster.sh inputFile
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/formatted/*; do qsub 06b_cluster.sh $i; done
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06b_cluster.sh "${fileList[$i]}"; done
## job 814062 to 814072
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_s4q20/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06b_cluster.sh "${fileList[$i]}"; done
## 869631 to 869641
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06b_cluster.sh "${fileList[$i]}"; done
## jobs 874194 to 874204

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
clusterOut=$outputsPath"/clustered_size_500_"$analysisTag
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
clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --cluster-size=500

# status message
echo "Analysis complete!"
