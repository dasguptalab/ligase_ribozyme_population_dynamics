#!/bin/bash

# script to run job scripts that count the number of sequences in sequence families
# usage: bash 13b_conservation.sh
## run 6
## jobs 1658215 to 1658222
## test_22May2025
## jobs 1728249 to 1728256

# retrieve analysis outputs absolute path
#outputsPath="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs"
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/test_22May2025/13b_overhang_conservation_above2"

# create outputs directory
mkdir $outDir

# read in sequence count data for the specified round
countsFile=$outputsPath"/09b_quantified_above2/counts_plot_table_noDoped.csv"

# loop over each input run num
for runNum in {1..8}; do 
	# retrieve the input round number
	roundNum=$runNum
	# status message
	echo "Beginning analysis of round $roundNum ..."
	# submit job script
	qsub 13_conserved.sh $roundNum $outDir $countsFile	
done

# status message
echo "Analysis complete!"
