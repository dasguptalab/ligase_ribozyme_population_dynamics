#!/bin/bash

# script to filter.fastq files and keep sequences with matching up- and down-stream sequences
# usage: bash 04_filter.sh 

# primer: GGCUAAGG -> GGCTAAGG
# library: GACUCACUGACACAGAUCCACUCACGGACAGCGG(Nx40)CGCUGUCCUUUUUUGGCUAAGG -> 96bp total
# target trimmed -> GGACAGCG(Nx40)CGCTGTCC(NxM) -> at least 56bp total

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
filterOut=$outputsPath"/filtered"
mkdir $filterOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $filterOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $filterOut

# status message
echo "Beginning analysis..."

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in $outputsPath"/trimmed/"*_uForward\.fq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/\.fq//')
	# status message
	echo "Processing $sampleTag"
	# filter to keep sequences with matching up- and down-stream sequences
	cat $f1 | grep -Ex -B1 '.*GGACAGCG.{40}CGCTGTCC.*' | sed "s/^.*GGACAGCG//g" | sed "s/CGCTGTCC.*$//g" | grep -Ex -B1 '.{40}' | grep -v "^--$" > $filterOut"/"$sampleTag".flt.fa"
	#cat $f1 | awk '/GGACAGCG.{40}CGCTGTCC/{if (a && a !~ /GGACAGCG.{40}CGCTGTCC/) print a; print} {a=$0}' | sed "s/^.*GGACAGCG//g" | sed "s/CGCTGTCC.*$//g" | grep -Ex -B1 '.{40}' > $filterOut"/"$sampleTag".flt.fa"
	# status message
	echo "$sampleTag processed!"
done

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in $outputsPath"/merged/"*\.fq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/\.fq//')
	# status message
	echo "Processing $sampleTag"
	# filter to keep sequences with matching up- and down-stream sequences
	cat $f1 | grep -Ex -B1 '.*GGACAGCG.{40}CGCTGTCC.*' | sed "s/^.*GGACAGCG//g" | sed "s/CGCTGTCC.*$//g" | grep -Ex -B1 '.{40}' | grep -v "^--$" > $filterOut"/"$sampleTag".flt.fa"
	#cat $f1 | awk '/GGACAGCG.{40}CGCTGTCC/{if (a && a !~ /GGACAGCG.{40}CGCTGTCC/) print a; print} {a=$0}' | sed "s/^.*GGACAGCG//g" | sed "s/CGCTGTCC.*$//g" | grep -Ex -B1 '.{40}' > $filterOut"/"$sampleTag".flt.fa"
	# status message
	echo "$sampleTag processed!"
done

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in $outputsPath"/merged/logs/"*_stiched_reads.fqiled\.fq_*\.fastq; do
	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/\.fastq//')
	# status message
	echo "Processing $sampleTag"
	# filter to keep sequences with matching up- and down-stream sequences
	cat $f1 | grep -Ex -B1 '.*GGACAGCG.{40}CGCTGTCC.*' | sed "s/^.*GGACAGCG//g" | sed "s/CGCTGTCC.*$//g" | grep -Ex -B1 '.{40}' | grep -v "^--$" > $filterOut"/"$sampleTag".flt.fa"
	#cat $f1 | awk '/GGACAGCG.{40}CGCTGTCC/{if (a && a !~ /GGACAGCG.{40}CGCTGTCC/) print a; print} {a=$0}' | sed "s/^.*GGACAGCG//g" | sed "s/CGCTGTCC.*$//g" | grep -Ex -B1 '.{40}' > $filterOut"/"$sampleTag".flt.fa"
	# status message
	echo "$sampleTag processed!"
done

# status message
echo "Analysis complete!"
