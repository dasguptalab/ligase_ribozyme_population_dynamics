#!/bin/bash

# script to filter sequences by length
# usage: bash 07_filterByLength.sh 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set the inputs directory
filterPairOut=$outputsPath"/filterByPair"

# set the ouputs directory
filterLengthOut=$outputsPath"/filterByLength"

# create the outputs directory
mkdir $filterLengthOut

# move to the new directory
cd $filterLengthOut

# loop over each file
for f1 in $filterPairOut; do
	# retrieve sample name
	sampleName=$(echo $f1 | sed "s/\.flt\.fq//g")
	# status message
	echo "Processing $sampleName ..."
	# filter to keep reads > 36 residues
	#cat $clusterOut"/combined.fmt.fasta" | sed "s/^>/>$(printf '%.0s~' {0..36})/g" | awk 'length($0) > 36' | sed "s/~//g" | grep -B1 --no-group-separator -v ">" > $clusterOut"/combined.flt.fmt.fasta"
	# filter to keep reads = 40 residues
	cat $f1 | sed "s/^>/>$(printf '%.0s~' {0..39})/g" | awk 'length($0) > 39' | sed "s/~//g" | awk 'length($0) < 41' | grep -B1 --no-group-separator -v ">" > $filterLengthOut"/"$sampleName".flt40.fasta"
done

# status message
echo "Analysis complete!"