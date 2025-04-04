#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_jobOutput
#$ -pe smp 8
#$ -q largemem

# script to cluster sequences using clustalo
# usage: qsub 07_cluster.sh sampleTag
# above 2
# 1500 -> default
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_above2\.fa; do echo $i; qsub 07_cluster.sh $i; done
## jobs 1541791 to 1541805
# 1400
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_above2\.fa; do echo $i; qsub 07_cluster.sh $i; done
## jobs 1545201 to 1545256
# all, 1500
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_formatted\.fa; do echo $i; qsub 07_cluster.sh $i; done
## jobs 1582635 to 1582645
# above 2, 1000
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_above2\.fa; do echo $i; qsub 07_cluster.sh $i; done
## jobs

# load the software module
module load bio/0724

# set the sample file
sampleFile=$1

# retrieve the sample tag
sampleTag=$(basename $sampleFile | sed 's/_formatted_above2\.fa//')
#sampleTag=$(basename $sampleFile | sed 's/_formatted\.fa//')
sampleTag=$sampleTag

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
outputsPath=$outputsPath"/07_clustered_1000"
mkdir $outputsPath

# move to the new directory
cd $outputsPath

# status message
echo "Processing $sampleFile ..."

# cluster sequences
#clustalo --threads=$NSLOTS -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=500 
#clustalo -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=500 --full --percent-id --distmat-out=$outputsPath"/"$sampleTag"_distances.txt"
clustalo --threads=$NSLOTS -i $inputsPath"/"$sampleFile --clustering-out=$outputsPath"/"$sampleTag"_clustered.aux" -o $outputsPath"/"$sampleTag"_aligned.fa" --cluster-size=1000

# status message
echo "Analysis complete!"
