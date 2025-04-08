#!/bin/bash

# script to run job scripts that count the number of sequences in sequence families
# usage: bash 11_identification.sh

# retrieve analysis outputs absolute path
#outputsPath="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs"
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
#outDir="/11a_family_identification"
outDir=$outputsPath"/11a_family_identification_above2"

# create outputs directory
mkdir $outDir

# read in cluster family sequence data
#peaksFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv"
peaksFile=$outputsPath"/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv"

# read in sequence count data for the specified round
#seqsFile="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09b_quantified_above2/counts_plot_table.csv"
seqsFile=$outputsPath"/09b_quantified_above2/counts_plot_table.csv"

# loop over each input run num
for runNum in {1..8}; do 
	# retrieve the input round number
	roundNum=$runNum
	# status message
	echo "Beginning analysis of round $roundNum ..."
	# submit job script
	qsub 11a_identification.sh $roundNum $outDir $peaksFile $seqsFile	
done

# status message
echo "Analysis complete!"
