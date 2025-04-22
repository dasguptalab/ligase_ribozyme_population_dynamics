#!/bin/bash

# script to run job scripts that count the number of sequences in sequence families
# usage: bash 11a_identification.sh
## run 6
## jobs 1621144 to 1621151

# retrieve analysis outputs absolute path
#outputsPath="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs"
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/11a_family_identification_all"

# create outputs directory
mkdir $outDir

# read in cluster family sequence data
peaksFile=$outputsPath"/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv"

# read in sequence count data for the specified round
seqsFile=$outputsPath"/09a_quantified_all/counts_plot_table.csv"
countsFile=$outputsPath"/09a_quantified_all/counts_plot_table_noDoped.csv"

# create the combined counts plotting data
head -1 $outputsPath"/09a_quantified_all/r1_S1_L001_in_r1_S1_L001_counts_plot_table.csv" > $seqsFile
for i in $outputsPath"/09a_quantified_all/"*_counts_plot_table.csv; do cat $i | tail -n+2 >> $seqsFile; done
cat $seqsFile | grep -v "doped" > $countsFile

# create the split sequence data files for each round


# loop over each input run num
for runNum in {1..8}; do 
	# retrieve the input round number
	roundNum=$runNum
	# status message
	echo "Beginning analysis of round $roundNum ..."
	# submit job script
	qsub 11_identify.sh $roundNum $outDir $peaksFile $countsFile	
done

# status message
echo "Analysis complete!"
