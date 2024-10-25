#!/bin/bash

# script to perform trimmomatic trimming of paired end reads
# usage: bash 02_trim.sh 

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# load the software module
module load bio/2.0

# set phred score for trimming
# https://www.drive5.com/usearch/manual/quality_score.html
score=33

# retrieve adapter absolute path for alignment
#adapterPath=$(grep "adapter:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/adapter://g")

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set inputs path
inputsPath=$outputsPath"/merged"

# make a new directory for analysis
outputsPath=$outputsPath"/trimmed"
mkdir $outputsPath
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $outputsPath

# status message
echo "Analyzing combined data..."

# unzip any gz read files
gunzip -v $inputsPath"/"*\.gz

# loop over un-filtered merged reads for each run
for f1 in $inputsPath"/"*\.extendedFrags.fastq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/\.extendedFrags\.fastq//')
	# status message
	echo "Processing $sampleTag"
	# perform adapter trimming on paired reads
	# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
	trimmomatic PE -threads 4 -phred"$score" $inputsPath"/"$sampleTag".notCombined_1.fastq" $inputsPath"/"$sampleTag".notCombined_2.fastq" $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" MAXINFO:56:0.8 ILLUMINACLIP:"$adapterPath":2:30:10:1:TRUE
	trimmomatic SE -threads 4 -phred"$score" $f1 $sampleTag"_trimmed.fq.gz" MAXINFO:56:0.8 ILLUMINACLIP:"$adapterPath":2:30:10:1:TRUE
	# status message
	echo "$sampleTag processed!"
done

# status message
echo "Analysis complete!"
