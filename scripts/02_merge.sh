#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_merge_NGmerge_jobOutput
#$ -pe smp 8
#$ -q largemem

# script to perform merging of paired end reads into single reads
# usage: qsub 02_merge.sh trimmed
## outputs_4quality15
## job 807481
## job 808607
## outputs_81quality20
## job 817815
## outputs_10quality20
## job 818009
# usage: qsub 02_merge.sh trimmed_s4q20
## job 817813
## job 
# usage: qsub 02_merge.sh trimmed_q20
## job 819645
## job

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# retrieve input analysis type
analysisType=$1

# retrieve analysis tag
analysisTag=$(echo $1 | sed "s/trimmed_//g")

# retrieve software absolute path
softwarePath=$(grep "software_NGmerge:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/software_NGmerge://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/"$analysisType

# make a new directory for analysis
mergeOut=$outputsPath"/merged_"$analysisTag
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
for f1 in $inputsPath"/"*_pForward\.fq*; do
	# set paired file name
	f2=$(echo $f1 | sed 's/_pForward/_pReverse/')
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_pForward\.fq.*//')
	# status message
	echo "Processing $sampleTag"
	#./NGmerge -v -n 8 -1 $f1 -2 $f2 -o $mergeOut"/"$sampleTag"_stiched_reads.fq" -m 8 -p 0 -d -e 8 -l $mergeOut"/logs/"$sampleTag"_log_stitching_results.txt" -f $mergeOut"/logs/"$sampleTag"_stiched_reads_failed.fq" -j $mergeOut"/logs/"$sampleTag"_log_formatted_alignments.txt" -q 33 -u 40
	./NGmerge -v -n 8 -1 $f1 -2 $f2 -o $mergeOut"/"$sampleTag"_stiched_reads.fq" -m 20 -p 0.1 -l $mergeOut"/logs/"$sampleTag"_log_stitching_results.txt" -f $mergeOut"/"$sampleTag"_stiched_reads_failed.fq" -j $mergeOut"/logs/"$sampleTag"_log_formatted_alignments.txt" -q 33 -u 40
	# status message
	echo "$sampleTag processed!"
done

#Print status message
echo "Analysis complete!"
