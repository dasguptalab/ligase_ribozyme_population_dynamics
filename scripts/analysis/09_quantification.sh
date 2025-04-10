#!/bin/bash

# script to run job scripts that count the number of sequences shared across runs
# usage: bash 09_quantification.sh
## quantification of sequences with at least 3 reads
## run 4
## jobs 1541857 to 1541922
## outputs
## jobs 1580433 to 1580579
## test_26March2025
## jobs 1581231 to 1581377
## run6
## jobs 1609679 to 1609838
## jobs 1614618 to 1614633


# loop over each input run num
for runNum in {1..8}; do 
	# loop over each run data num
	for dataNum in {1..8}; do 
		# setup run name tags
		runInput="r"$runNum"_S"$runNum"_L001"
		runData="r"$dataNum"_S"$dataNum"_L001"
		# status message
		echo "Beginning analysis of $runInput over $runData ..."
		# submit job script
		qsub 09a_quantify.sh $runInput $runData
		qsub 09b_quantify.sh $runInput $runData
	done
	# status message
	echo "Beginning analysis of $runInput top 10 sequences..."
	# submit job scripts for the top 10 sequences per round
	qsub 09c_quantify.sh $runInput
	qsub 09d_quantify.sh $runInput
done

# combine run data
#head -1 /scratch365/ebrooks5/RNA_evolution/outputs/09a_quantified/r8_S8_L001_in_r8_S8_L001_counts_plot_table.csv > /scratch365/ebrooks5/RNA_evolution/outputs/09a_quantified/counts_plot_table.csv
#for i in /scratch365/ebrooks5/RNA_evolution/outputs/09a_quantified/*_counts_plot_table.csv; do echo $i; tail -n+2 $i >> /scratch365/ebrooks5/RNA_evolution/outputs/09a_quantified/counts_plot_table.csv; done
#cat /scratch365/ebrooks5/RNA_evolution/outputs/09a_quantified/counts_plot_table.csv | grep -v "doped" > /scratch365/ebrooks5/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv

# status message
echo "Analysis complete!"
