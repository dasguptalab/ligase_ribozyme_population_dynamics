#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N cluster_families_jobOutput
#$ -q largemem

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: qsub 10_cluster_families.sh inputFile
# usage ex: qsub 10_cluster_families.sh 07a_clustered r8_S8_L001
## job 985948
# usage ex: qsub 10_cluster_families.sh 07a_clustered r7_S7_L001
## job 
# usage ex: qsub 10_cluster_families.sh 07a_clustered r6_S6_L001
## job 
# usage ex: qsub 10_cluster_families.sh 07a_clustered doped21-r3_S12_L001
## job
# usage ex: qsub 10_cluster_families.sh 07a_clustered doped21-r2_S11_L001
## job

# retrieve input file
inputFile=$1

# retrieve input sequences
inputSeqs=$2"_formatted_above9_cluster_sequences_identity_table.csv"

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/09_cluster_identity/"$inputFile

# make a new directory for analysis
tablesOut=$outputsPath"/10_cluster_families"
mkdir $tablesOut

# make a new directory for analysis
tablesOut=$tablesOut"/"$inputFile
mkdir $tablesOut
# check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $tablesOut directory already exsists... please remove before proceeding."
#	exit 1
#fi

# move to outputs directory
cd $tablesOut

# name output file
countsOut=$tablesOut"/"$2"_cluster_families_table.csv"

# retrieve header
header=$(head -1 $inputSeqs)

# add the run names to the header 
header=$(echo $header",doped21-r1_counts,doped21-r2_counts,doped21-r3_counts,r1_counts,r2_counts,r3_counts,r4_counts,r5_counts,r6_counts,r7_counts,r8_counts")

# add a header to the counts data tmp file
echo $header > $countsOut".tmp.csv"

# status message
echo "Beginning analysis..."
	
# loop over cluster sequences
while read seqData; do
	# retrieve the seq
	seq=$(echo $seqData | cut -d"," -f5)
	# reverse complement
	seq=$(echo $seq | tr ACGTacgt TGCAtgca | rev)
	# update the count data
	countData=$seqData
	# status message
	echo "Processing $seq ..."
	# loop over each round sequences file
	for f2 in $outputsPath"/05_combined/"*_combined\.fa; do
		# retrieve run name
		#runName=$(basename $f2 | cut -d"_" -f1)
		# status message
		#echo "Processing $runName ..."
		# count the number of seq occurances in each round
		numReads=$(cat $f2 | grep -wc $seq)
		# add the number of seqs for the round
		countData=$(echo $countData","$numReads)
	done
	# add the counts data to the tmp file
	echo $countData >> $countsOut".tmp.csv"
done < <(tail -n+2 $inputsPath"/"$inputSeqs)

# output the header and run counts
echo $header > $countsOut

# output the counts data
tail -n+2 $countsOut".tmp.csv" >> $countsOut

# clean up
rm $countsOut".tmp.csv"

# status message
echo "Analysis complete!"
