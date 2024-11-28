#!/bin/bash

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: bash 08_summarize.sh inputFile
# usage ex: bash 08_summarize.sh 07a_clustered
# usage ex: bash 08_summarize.sh 07b_clustered

# retrieve input file
inputFile=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/"$inputFile

# make a new directory for analysis
tablesOut=$outputsPath"/08_summarized"
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
for f1 in $inputsPath"/"*_clustered\.aux; do
	# trim to sample tag
	nameTag=$(basename $f1 | sed 's/_clustered\.aux//')
	# status message
	echo "Processing $nameTag ..."
	# retrieve cluster names
	cat $f1 | cut -d":" -f1 | sed "s/ /_/g" > $tablesOut"/"$nameTag"_cluster_numbers.tmp.txt"
	# retrieve sequence headers
	cat $f1 | cut -d":" -f2 | tr -s " " | cut -d" " -f8 > $tablesOut"/"$nameTag"_seq_headers.tmp.txt"
	# combine cluster numbers with sequence headers
	paste -d "," $tablesOut"/"$nameTag"_seq_headers.tmp.txt" $tablesOut"/"$nameTag"_cluster_numbers.tmp.txt" > $tablesOut"/"$nameTag"_cluster_sequences_table.tmp.csv"
	# add header to the output sequences file
	echo "run_name,sequence_ID,read_counts,cluster_ID,sequence" > $tablesOut"/"$nameTag"_cluster_sequences_table.clust.tmp.csv"
	# loop over each sequence line and retrieve sequence data
	while read line; do
		# retrieve sequence header
		seqHeader=$(echo $line | cut -d"," -f1-3)
		# retrieve sequence data
		# reverse complement seqeuence -> tr ACGTacgt TGCAtgca | rev
		seqData=$(cat $outputsPath"/06_formatted/"$nameTag".fa" | grep -A1 $seqHeader | tail -1 | tr ACGTacgt TGCAtgca | rev)
		# add sequence data to the output file
		echo $line","$seqData >> $tablesOut"/"$nameTag"_cluster_sequences_table.clust.tmp.csv"
	done < $tablesOut"/"$nameTag"_cluster_sequences_table.tmp.csv"
	# filter to keep only the peak sequences
	cat $tablesOut"/"$nameTag"_cluster_sequences_table.clust.tmp.csv" | awk -F',' '!seen[$4]++' > $tablesOut"/"$nameTag"_cluster_peaks_table.tmp.csv"
	# remove the header
	cat $tablesOut"/"$nameTag"_cluster_peaks_table.tmp.csv" | tail -n+2 > $tablesOut"/"$nameTag"_cluster_peaks_table.noHeader.tmp.csv"
	# add header to the output sequences file
	echo "run_name,sequence_ID,read_counts,cluster_ID,sequence_counts,sequence" > $tablesOut"/"$nameTag"_cluster_peaks_table.tmp.csv"
	# loop over each peak line and add the number of sequences in the cluster
	while read line; do
		# retrieve cluster number
		clusterName=$(echo $line | cut -d"," -f4)
		# count the number of sequences in the current cluster
		clusterCount=$(cat $tablesOut"/"$nameTag"_cluster_sequences_table.clust.tmp.csv" | grep -c $clusterName)
		# separate sequence from the other data
		seqHeader=$(echo $line | cut -d"," -f1-4)
		seqData=$(echo $line | cut -d"," -f5)
		# add cluster count to outputs
		echo $seqHeader","$clusterCount","$seqData >> $tablesOut"/"$nameTag"_cluster_peaks_table.tmp.csv"
	done < $tablesOut"/"$nameTag"_cluster_peaks_table.noHeader.tmp.csv"
	# re-format cluster names
	cat $tablesOut"/"$nameTag"_cluster_sequences_table.clust.tmp.csv" | sed "s/Cluster_//g" > $tablesOut"/"$nameTag"_cluster_sequences_table.csv"
	cat $tablesOut"/"$nameTag"_cluster_peaks_table.tmp.csv" | sed "s/Cluster_//g" > $tablesOut"/"$nameTag"_cluster_peaks_table.csv"
	# clean up
	rm $tablesOut"/"$nameTag"_cluster_numbers.tmp.txt"
	rm $tablesOut"/"$nameTag"_seq_headers.tmp.txt"
	rm $tablesOut"/"$nameTag"_cluster_sequences_table"*".tmp.csv"
	rm $tablesOut"/"$nameTag"_cluster_peaks_table"*".tmp.csv"
done

# combine tables from each run
cat $tablesOut"/"*"_cluster_sequences_table.csv" | head -1 > $tablesOut"/cluster_sequences_table.csv"
for i in $tablesOut"/"*"_cluster_sequences_table.csv"; do tail -n+2 $i >> $tablesOut"/cluster_sequences_table.csv"; done
cat $tablesOut"/"*"_cluster_peaks_table.csv" | head -1 > $tablesOut"/cluster_peaks_table.csv"
for i in $tablesOut"/"*"cluster_peaks_table.csv"; do tail -n+2 $i >> $tablesOut"/cluster_peaks_table.csv"; done

# status message
echo "Analysis complete!"
