#!/bin/bash

# script to filter sequences by length
# usage: bash 07_filterByLength.sh 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set directory for inputs
formatOut=$outputsPath"/formatted"

cat $clusterOut"/combined.fasta" | sed "s/^>/>$(printf '%.0s~' {0..37})/g" | awk 'length($0) > 36' | sed "s/~//g" > $clusterOut"/combined.flt.fasta"
