#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_merge_NGmerge_jobOutput
#$ -pe smp 8

# script to perform merging of paired end reads into single reads
# usage: qsub 03_merge.sh

# retrieve input analysis type
analysisType=$1

# retrieve software absolute path
softwarePath=$(grep "software_NGmerge:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/software_NGmerge://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve trimmed reads path
readPath=$outputsPath"/trimmed"

# make a new directory for analysis
mergeOut=$outputsPath"/merged"
mkdir $mergeOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $mergeOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the software directory
cd $softwarePath

# status message
echo "Processing..."

# loop through all forward and reverse reads and merge each pair into a single read
for f1 in $readPath"/"*_pForward\.fq\.gz; do
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_pForward\.fq\.gz//')
	# set paired file name
	f2=$curSample"_pReverse.fq.gz"
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_pForward\.fq\.gz//')
	# status message
	echo "Processing $sampleTag"
	./NGmerge -v -n 8 -1 $f1 -2 $f2 -o $mergeOut"/stiched_reads.fa" -m 20 -p 0 -l $mergeOut"/log_stitching_results.txt" -f $mergeOut"/stiched_reads_failed.fa" -j $mergeOut"/log_formatted_alignments.txt" -q 33 -u 40
	# status message
	echo "$sampleTag processed!"
done

#Print status message
echo "Analysis complete!"
