#!/bin/bash

# script to perform trimmomatic trimming of paired end reads
# usage: bash 01a_trim.sh 

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# load the software module
module load bio/2.0

# retrieve paired reads absolute path for alignment
inputsPath=$(grep "pairedReads:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/pairedReads://g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve adapter absolute path for alignment
#adapterPath=$(grep "adapter:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/adapter://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
outputsPath=$outputsPath"/trimmed_"$analysisTag
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsPath

# set phred score for trimming
# https://www.drive5.com/usearch/manual/quality_score.html
score=33

# status message
echo "Beginning analysis..."

# loop through all samples
for f1 in $inputsPath"/"*_R1_001\.fastq\.gz; do
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R1_001\.fastq\.gz//')
	# set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_R1_001\.fastq\.gz//')
	# status message
	echo "Processing $sampleTag"
	# perform adapter trimming on paired reads
	# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
	#trimmomatic PE -threads 4 -phred"$score" $f1 $f2 $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath":2:30:10 SLIDINGWINDOW:4:15 HEADCROP:23 MINLEN:56
	#trimmomatic PE -threads 4 -phred"$score" $f1 $f2 $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" AVGQUAL:20
	trimmomatic PE -threads 4 -phred"$score" $f1 $f2 $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath":2:30:10:1:TRUE MAXINFO:40:0.8
	# status message
	echo "$sampleTag processed!"
done

# status message
echo "Analysis complete!"
