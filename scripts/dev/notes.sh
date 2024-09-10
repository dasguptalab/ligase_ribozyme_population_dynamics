#!/bin/bash

# set inDir
inDir="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/filteredByRegionLength"

# set outDir
outDir="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/analyzedByRegionLength_test"

for i in $inDir"/"*"flt.fq"; do 
	# retrieve sample name
	sampleName=$(basename $i | sed "s/\flt.\fq//g")
	# sattus message
	echo "Processing $sampleName ..."
	# generate stats
	cat $i | awk 'NR % 2 == 0' | sort -n | uniq -c | sort -rk1 > $outDir"/"$sampleName"_seqStats.flt.txt"
	# count number of unique sequences
	unique=$(wc -l $outDir"/"$sampleName"_seqStats.flt.txt")
	# print number of unique sequences
	echo $unique
done
