#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_clustalo_jobOutput
#$ -pe smp 8
#$ -q largemem

# script to cluster sequences using clustalo
# usage: qsub 06a_cluster.sh inputFile
# usage ex: qsub 06a_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs/formatted/r8_S8_L001_combined.fmt.fa
## job 816015 -> FATAL: Memory allocation for distance matrix failed
# usage ex: qsub 06a_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs/formatted/r7_S7_L001_combined.fmt.fa
## job 816016 -> FATAL: Memory allocation for distance matrix failed
# usage ex: qsub 06a_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs/formatted/r6_S6_L001_combined.fmt.fa
## job 816017
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/formatted/*; do qsub 06a_cluster.sh $i; done
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06a_cluster.sh "${fileList[$i]}"; done
## job 814037 to 814047
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_s4q20/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06a_cluster.sh "${fileList[$i]}"; done
## job 868912 to 868923
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06a_cluster.sh "${fileList[$i]}"; done
## jobs 874029 to 874039
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 06a_cluster.sh "${fileList[$i]}"; done
## jobs 880073 to 880083

# load the software module
module load bio/0724

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# clean up input file name
nameTag=$(basename $inputFile | sed "s/\./_/g" | sed "s/_fa/\.fa/g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# make a directory for the clustering outputs
clusterOut=$outputsPath"/clustered_"$analysisTag
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
clustalo --threads=$NSLOTS -i $inputFile --clustering-out=$clusterOut"/"$nameTag"_clustered.aux" -o $clusterOut"/"$nameTag"_aligned.fa"

# status message
echo "Analysis complete!"
