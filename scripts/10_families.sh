#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_families_jobOutput
#$ -q largemem

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: qsub 10_families.sh inputFile inputRun
# usage ex: qsub 10_families.sh 07a_clustered ...
# r8_S8 to doped21-r2_S11
## job 985948 to 985999
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/08_summarized/*_formatted_above9_cluster_sequences_table.csv; do runInput=$(basename $i | sed "s/_formatted_above9_cluster_sequences_table\.csv//g"); qsub 10_families.sh 07a_clustered $runInput; done
# job 1007705 to 1007715

# retrieve input file
inputFile=$1

# retrieve input run
inputRun=$2

# retrieve input sequences
inputSeqs=$inputRun"_formatted_above9_cluster_sequences_table.csv"

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/08_summarized/"$inputFile

# make a new directory for analysis
tablesOut=$outputsPath"/10_families"
mkdir $tablesOut

# make a new directory for analysis
tablesOut=$tablesOut"/"$inputFile
mkdir $tablesOut

# move to outputs directory
cd $tablesOut

# name output file
countsOut=$tablesOut"/"$inputRun"_counts_table.csv"
countsPlotOut=$tablesOut"/"$inputRun"_counts_plot_table.csv"

# retrieve header
inputHeader=$(head -1 $inputsPath"/"$inputSeqs)

# add the run names to the header 
header=$(echo $inputHeader",doped21-r1_counts,doped21-r2_counts,doped21-r3_counts,r1_counts,r2_counts,r3_counts,r4_counts,r5_counts,r6_counts,r7_counts,r8_counts")
headerPlot=$(echo $inputHeader",counts,counts_run_name")

# add a header to the counts data tmp file
echo $header > $countsOut
echo $headerPlot > $countsPlotOut

# status message
echo "Beginning analysis of $inputRun ..."

# loop over round sequences
while read data; do
	# clean up the run name
	seqData=$(echo $data | sed "s/>r//g")
	# remove the sequence
	dataNoSeq=$(echo $seqData | cut -d"," -f1-4)
	# retrieve the seq
	seq=$(echo $seqData | cut -d"," -f5)
	# reverse compliment the sequence
	seqRev=$(echo $seq | tr ACGTacgt TGCAtgca | rev)
	# update the count data
	countData=$dataNoSeq","$seqRev
	# status message
	echo "Processing $seq ..."
	# loop over each round sequences file
	for f2 in $outputsPath"/05_combined/"*_combined\.fa; do
		# retrieve run name
		runName=$(basename $f2 | sed "s/_S.*_L001_combined\.fa//g" | sed "s/r//g" | sed "s/21-/_/g")
		# count the number of seq occurances in each round
		numReads=$(cat $f2 | grep -wc $seq)
		# add the number of seqs for the round
		countData=$(echo $countData","$numReads)
		# add the counts data to the tmp file
		echo $seqData","$numReads","$runName >> $countsPlotOut
	done
	# add the counts data to the tmp file
	echo $countData >> $countsOut
done < <(tail -n+2 $inputsPath"/"$inputSeqs)
	
# loop over cluster sequences
#while read seqData; do
	# retrieve the seq
#	seq=$(echo $seqData | cut -d"," -f5)
	# reverse complement
	#seq=$(echo $seq | tr ACGTacgt TGCAtgca | rev)
	# update the count data
#	countData=$seqData
	# status message
#	echo "Processing $seq ..."
	# loop over each round sequences file
#	for f2 in $outputsPath"/05_combined/"*_combined\.fa; do
	#for f2 in $outputsPath"/06_formatted/"*_formatted\.fa; do
		# count the number of seq occurances in each round
#		numReads=$(cat $f2 | grep -wc $seq)
		# add the number of seqs for the round
#		countData=$(echo $countData","$numReads)
#	done
	# add the counts data to the tmp file
#	echo $countData >> $countsOut".tmp.csv"
#done < <(tail -n+2 $inputsPath"/"$inputSeqs)

# output the header and run counts
#echo $header > $countsOut

# output the counts data
#tail -n+2 $countsOut".tmp.csv" >> $countsOut

# clean up
#rm $countsOut".tmp.csv"

# status message
echo "Analysis complete!"
