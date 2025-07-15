#!/bin/bash

# script to run job scripts that count the number of sequences in sequence families
# usage: bash 13a_conservation.sh

# retrieve analysis outputs absolute path
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/13d_overhang_conservation_families"

# create outputs directory
mkdir $outDir

# read in sequence count data for the specified round
countsFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/ST2_family_table/r8_family_count_data.csv"

# run the analysis
Rscript 13_sequence_overhang_conservation.R $outDir $seqsFile

# status message
echo "Analysis complete!"
