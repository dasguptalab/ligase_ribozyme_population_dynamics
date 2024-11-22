#!/bin/bash

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: bash 09_cluster_identity.sh inputFile
# usage ex: bash 09_cluster_identity.sh 07a_clustered
# usage ex: bash 09_cluster_identity.sh 07b_clustered

# retrieve input file
inputFile=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/08_cluster_tables/"$inputFile

# make a new directory for analysis
tablesOut=$outputsPath"/09_cluster_identity"
mkdir $tablesOut

# make a new directory for analysis
tablesOut=$tablesOut"/"$inputFile
mkdir $tablesOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $tablesOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to outputs directory
cd $tablesOut

# status message
echo "Beginning analysis..."

# loop through all samples
for f1 in $inputsPath"/"*_above9_cluster_peaks_table\.csv; do
	# trim to sample tag
	nameTag=$(basename $f1 | sed 's/_peaks_table\.csv//g')
	# setup peak sequences file name
	f2=$inputsPath"/"$nameTag"_sequences_table.csv"
	# status message
	echo "Processing $nameTag ..."
	# initialize counter
	peakCount=0
	# loop over each cluster name number
	while read peakData; do
		# increment peak counter
		peakCount=$(($peakCount+1))
		# check if header line
		if [ $peakCount -eq 1 ]; then
			# adjust and output the header
			echo $peakData",avg_ID,SD_ID,highest_ID,lowest_ID" > $tablesOut"/"$nameTag"_peaks_identity_table.csv"
			echo $peakData",percent_ID" > $tablesOut"/"$nameTag"_sequences_identity_table.csv"
		else
			# get the current cluster number
			clusterName=$(echo $peakData | cut -d"," -f4)
			# status message
			echo "Processing Cluster $clusterName ..."
			# retrieve sequences for the current cluster
			cat $f2 | grep ",$clusterName," > $tablesOut"/"$nameTag"_cluster_"$clusterName"_table.csv"
			# initialize counters
			numSeqs=0
			totalID=0
			totalSqrVar=0
			# loop over each sequence line and retrieve sequences
			while read seqData; do
				# retrieve the current cluster peak sequence
				peakSeq=$(echo $peakData | cut -d"," -f6)
				# retrieve the current sequence in the current cluster
				currSeq=$(echo $seqData | cut -d"," -f5)
				# initialize mismatch counter
				numMismatch=0
				# loop over each character of the sequence
				#echo $sequence | awk '{for (i=0; ++i <= length($0);) printf "%s", substr($0, i, 1)}'
				for (( i=0; i<${#currSeq}; i++ )); do
					# compare each character with the peak
					if [ "${currSeq:$i:1}" != "${peakSeq:$i:1}" ]; then
						# increment the number of mismatches
						numMismatch=$(($numMismatch+1))
					fi
				done
				# calculate the number of matches
				seqLength=40
				# calculate the percent identity
				distID=$(($seqLength - $numMismatch))
				propID=$(echo "scale=4; $distID / $seqLength" | bc)
				percentID=$(echo "scale=4; $propID * 100" | bc)
				# add the percent identity to the running total for the cluster
				totalID=$(echo "scale=4; $totalID + $percentID" | bc)
				# increase the counter for the number of sequences in the cluster
				numSeqs=$(($numSeqs + 1))
				# add the percent identity to the current sequence information and output to a new table
				echo $seqData","$percentID >> $tablesOut"/"$nameTag"_sequences_identity_table.csv"
			done < $tablesOut"/"$nameTag"_cluster_"$clusterName"_table.csv"
			# calculate the avg percent identity for the cluster
			avgID=$(echo "scale=4; $totalID / $numSeqs" | bc)
			# reset counter
			lineCount=0
			# loop over each sequence line and calculate the sd percent identity for the cluster
			while read seqID; do
				# increment line counter
				lineCount=$(($lineCount + 1))
				# check if not the header line
				if [ $lineCount -ne 1 ]; then
					# retrieve the current percent identity
					currID=$(echo $seqID | cut -d"," -f7)
					# begin to calculate the sd
					varID=$(echo "scale=4; $currID - $avgID" | bc)
					sqrVar=$(echo "scale=4; $varID * $varID" | bc)
					totalSqrVar=$(echo "scale=4; $totalSqrVar + $sqrVar" | bc)
				fi
			done < $tablesOut"/"$nameTag"_sequences_identity_table.csv"
			# complete the sd calculation
			divSD=$(echo "scale=4; $totalSqrVar / $numSeqs" | bc)
			sdID=$(bc <<< "scale=4; sqrt($divSD)")
			# find the highest and lowest values
			highest=$(cat $tablesOut"/"$nameTag"_sequences_identity_table.csv" | cut -d"," -f6 | sort -rn | head -1)
			lowest=$(cat $tablesOut"/"$nameTag"_sequences_identity_table.csv" | cut -d"," -f6 | sort -rn | sed '$d' | tail -1)
			# add the avg and SD percent identity to the peak sequence information and output to a new table
			echo $peakData","$avgID","$sdID","$highest","$lowest >> $tablesOut"/"$nameTag"_peaks_identity_table.csv"
			# clean up
			rm $tablesOut"/"$nameTag"_cluster_"$clusterName"_table.csv"
			# sanity checks
			echo "avgID: "$avgID
			echo "sdID: "$sdID
			echo "highest: "$highest
			echo "lowest: "$lowest
		fi
	done < $f1
done

# combine tables from each run
cat $tablesOut"/"*"_sequences_identity_table.csv" | head -1 > $tablesOut"/sequences_identity_table.csv"
for i in $tablesOut"/"*"_sequences_identity_table.csv"; do tail -n+2 $i >> $tablesOut"/sequences_identity_table.csv"; done
cat $tablesOut"/"*"_peaks_identity_table.csv" | head -1 > $tablesOut"/peaks_identity_table.csv"
for i in $tablesOut"/"*"_peaks_identity_table.csv"; do tail -n+2 $i >> $tablesOut"/peaks_identity_table.csv"; done

# status message
echo "Analysis complete!"
