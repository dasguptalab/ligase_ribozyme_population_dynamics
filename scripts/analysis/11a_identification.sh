#!/bin/bash

# script to run job scripts that count the number of sequences in sequence families
# usage: bash 11a_identification.sh
## run 6
## jobs 1621144 to 1621151
## tests
### run 1
## jobs 1689673 to 1689787

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
seqsFileAll=$outputsPath"/09a_quantified_all/counts_plot_table.csv"
seqsFile=$outputsPath"/09a_quantified_all/counts_plot_table_noDoped.csv"

# status message
#echo "Preparing data for analysis..."

# create the combined counts plotting data
#head -1 $outputsPath"/09a_quantified_all/r1_S1_L001_in_r1_S1_L001_counts_plot_table.csv" > $seqsFileAll
#for i in $outputsPath"/09a_quantified_all/"*_counts_plot_table.csv; do cat $i | tail -n+2 >> $seqsFileAll; done
#cat $seqsFileAll | grep -v "doped" > $seqsFile

# status message
#echo "Data prepared!"

# loop over each input run num
#for runNum in {1..8}; do 
	# retrieve the input round number
	#roundNum=$runNum
	roundNum=1
	# status message
	echo "Beginning analysis of round $roundNum ..."
	# create run name tage
	runTag="r"$roundNum"_S"$roundNum"_L001"
	# subset the seqs file for the current round
	roundFile=$outDir"/"$runTag"_counts_plot_table_noDoped.tmp.csv"
	cat $seqsFile | grep $runTag > $roundFile
	# create the split sequence data files for each round
	inputSeqsData=$outDir"/"$runTag"_counts_plot_table_noDoped.tmp"
	split --lines=50000 $roundFile $inputSeqsData".split.fa"
	# clean up
	rm $roundFile
	# initialize split counter
	splitNum=0
	# loop over each split file
	for inputSeqs in $inputSeqsData".split.fa"*; do
		# increment split counter
		splitNum=$(($splitNum+1))
		# submit job script
		qsub 11a_identify.sh $roundNum $outDir $peaksFile $inputSeqs $splitNum
	done
#done

# status message
echo "Analysis complete!"
