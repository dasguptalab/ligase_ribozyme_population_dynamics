#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_merge_NGmerge_jobOutput
#$ -pe smp 8
#$ -q largemem

# script to perform merging of paired end reads into single reads
# usage: qsub 03_merge.sh
## outputs_4quality15
## job 807481
## job 808607
## outputs_4quality20

## outputs_81quality20


# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# retrieve input analysis type
analysisType=$1

# retrieve software absolute path
softwarePath=$(grep "software_NGmerge:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/software_NGmerge://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/trimmed"

# make a new directory for analysis
mergeOut=$outputsPath"/merged"
mkdir $mergeOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $mergeOut directory already exsists... please remove before proceeding."
	exit 1
fi

# create output logs directory
mkdir $mergeOut"/logs"

# move to the software directory
cd $softwarePath

# status message
echo "Processing..."

# loop through all forward and reverse reads and merge each pair into a single read
for f1 in $inputsPath"/"*_pForward\.fq\.gz; do
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_pForward\.fq\.gz//')
	# set paired file name
	f2=$curSample"_pReverse.fq.gz"
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_pForward\.fq\.gz//')
	# status message
	echo "Processing $sampleTag"
	./NGmerge -v -n 8 -1 $f1 -2 $f2 -o $mergeOut"/"$sampleTag"_stiched_reads.fq" -m 8 -p 0 -d -e 8 -l $mergeOut"/logs/"$sampleTag"_log_stitching_results.txt" -f $mergeOut"/logs/"$sampleTag"_stiched_reads_failed.fq" -j $mergeOut"/logs/"$sampleTag"_log_formatted_alignments.txt" -q 33 -u 40
	# status message
	echo "$sampleTag processed!"
done

# unzip all paired reads
gunzip -v $mergeOut"/"*\.fq\.gz
gunzip -v $mergeOut"/logs/"*\.fastq\.gz

#Print status message
echo "Analysis complete!"
