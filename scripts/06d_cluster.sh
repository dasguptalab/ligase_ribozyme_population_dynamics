#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_cd_hit_jobOutput
#$ -pe smp 63
#$ -q largemem

# script to cluster sequences using cd-hit
# usage: qsub 06d_cluster.sh inputFile
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06d_cluster.sh "${fileList[$i]}"; done
## jobs 874882 to 874893
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/cleaned_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06d_cluster.sh "${fileList[$i]}"; done
## jobs 874898 to 874908
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06d_cluster.sh "${fileList[$i]}"; done
## jobs 874965 to 874975
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/cleaned_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06d_cluster.sh "${fileList[$i]}"; done
## jobs 874953 to 874963

# retrieve input file
inputFile=$1

# retrieve software path
softwarePath=$(grep "software_cdhit:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/software_cdhit://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# clean up input file name
nameTag=$(basename $inputFile | sed "s/\.fa//g" | sed "s/\./_/g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# make a directory for the clustering outputs
clusterOut=$outputsPath"/clustered_cd_hit_est_"$analysisTag
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
#cd-hit -T $NSLOTS -M 2000 -sc 1 -sf 1 -i $inputFile -o $clusterOut"/clustered_"$nameTag
$softwarePath"/"cd-hit-est -T $NSLOTS -M 2000 -c 0.95 -n 8 -i $inputFile -o $clusterOut"/clustered_"$nameTag

# status message
echo "Analysis complete!"
