#!/bin/bash

# script to perform bbduk trimming of paired end reads
# usage: bash 01b_trim.sh 

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# retrieve software path
softwarePath=$(grep "bbmap:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/bbmap://g")

# retrieve paired reads absolute path for alignment
inputsPath=$(grep "pairedReads:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/pairedReads://g")
# retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/adapter://g")
# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
trimOut=$outputsPath"/trimmed_q10_ftm5"
mkdir $trimOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $trimOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $softwarePath

# status message
echo "Beginning analysis..."

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in $inputsPath"/"*_R1_001\.fastq\.gz; do
	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R1_001\.fastq\.gz//')
	# set paired file name
	f2=$curSample"_R2_001.fastq.gz"
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_R1_001\.fastq\.gz//')
	# status message
	echo "Processing $sampleTag"
	# perform adapter trimming and quality filtering of paired reads
	./bbduk.sh qtrim=rl trimq=10 ftm=5 in=$f1 in2=$f2 out=$trimOut"/"$sampleTag"_pForward.fq" out2=$trimOut"/"$sampleTag"_pReverse.fq" outm=$trimOut"/"$sampleTag"_R1_001_failed.fq" outm2=$trimOut"/"$sampleTag"_R2_001_failed.fq" outs=$trimOut"/"$sampleTag"_unpaired.fq"
	# status message
	echo "$sampleTag processed!"
done

# status message
echo "Analysis complete!"
