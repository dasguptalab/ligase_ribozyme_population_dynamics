#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_jobOutput
#$ -pe smp 8

# script to cluster sequences using clustalo
# usage: qsub 07_cluster.sh sampleTag
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_formatted\.fa; do sampleTag=$(basename $i | sed 's/_formatted\.fa//'); echo $sampleTag; qsub 07_cluster.sh $sampleTag; done
## jobs 1267159 to 1267169
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted_above4/*_above4\.fa; do sampleTag=$(basename $i | sed 's/_formatted_above4\.fa//'); sampleTag=$sampleTag"_above4"; echo $sampleTag; qsub 07_cluster.sh $sampleTag $i; done
## jobs 1267177 to 1267187
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted_above2/*_above2\.fa; do sampleTag=$(basename $i | sed 's/_formatted_above2\.fa//'); sampleTag=$sampleTag"_above2"; echo $sampleTag; qsub 07_cluster.sh $sampleTag $i; done
## jobs 1502980 to 1502990
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted_above2/*_above2\.fa; do sampleTag=$(basename $i | sed 's/_formatted_above2\.fa//'); sampleTag=$sampleTag"_above2"; echo $sampleTag; qsub 07_cluster.sh $sampleTag $i; done
## jobs 1538063 to 1538074
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted_above2/*_above2\.fa; do sampleTag=$(basename $i | sed 's/_formatted_above2\.fa//'); sampleTag=$sampleTag"_above2"; echo $sampleTag; qsub 07_cluster.sh $sampleTag $i; done
## jobs

# load the software module
module load bio/0724

# retrieve the same tag
sampleTag=$1

# set the sample file
sampleFile=$2

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
outputsPath=$outputsPath"/07_clustered_1100_above2"
mkdir $outputsPath

# move to the new directory
cd $outputsPath

# status message
echo "Processing $sampleTag ..."

# cluster sequences
#clustalo --threads=$NSLOTS -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=500 
#clustalo -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=500 --full --percent-id --distmat-out=$outputsPath"/"$sampleTag"_distances.txt"
clustalo -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=1100

# status message
echo "Analysis complete!"
