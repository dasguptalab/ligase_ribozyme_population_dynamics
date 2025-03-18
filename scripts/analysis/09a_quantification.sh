#!/bin/bash

# script to run job scripts that count the number of sequences shared across runs
# usage: bash 09a_quantification.sh

# loop over each input run num
for runNum in {1..8}; do 
	# loop over each run data num
	for dataNum in {1..8}; do 
		# status message
		echo "Beginning analysis of $runNum over $dataNum ..."
		# setup run name tags
		runInput="r"$runNum"_S"$runNum"_L001"
		runData="r"$dataNum"_S"$dataNum"_L001"
		# submit job script
		qsub 09a_quantify.sh $runInput $runData
	done
done

# status message
echo "Analysis complete!"
