#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_quantify_jobOutput
#$ -q largemem

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: qsub 11_quantify.sh inputFile
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_formatted.fa; do runInput=$(basename $i | sed "s/_formatted.fa//g"); echo qsub 11_quantify.sh 07a_clustered $runInput; done
## job 
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_formatted.fa; do runInput=$(basename $i | sed "s/_formatted.fa//g"); echo qsub 11_quantify.sh 07b_clustered $runInput; done
## job 

# retrieve input file
inputFile=$1

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/06_formatted"

# retrieve input sequences
inputSeqs=$inputsPath"/"$2"_formatted.fa"

# make a new directory for analysis
tablesOut=$outputsPath"/11_quantified"
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

# name formatted sequences file
fmtSeqs=$inputsPath"/"$2"_formatted.tmp.fa"

# re-format input sequences for processing
cat $inputSeqs | tr "\n" "," | sed "s/>/\n>/g" | sed "s/,$//g" | sed '/^[[:space:]]*$/d' > $fmtSeqs
echo "" >> $fmtSeqs

# name output file
countsOut=$tablesOut"/"$2"_counts_table.csv"
countsPlotOut=$tablesOut"/"$2"_counts_plot_table.csv"

# retrieve header
inputHeader="run_name,sequence_ID,read_counts,sequence"

# add the run names to the header 
header=$(echo $inputHeader",doped21-r1_counts,doped21-r2_counts,doped21-r3_counts,r1_counts,r2_counts,r3_counts,r4_counts,r5_counts,r6_counts,r7_counts,r8_counts")
headerPlot=$(echo $inputHeader",counts,counts_run_name")

# add a header to the counts data tmp file
echo $header > $countsOut
echo $headerPlot > $countsPlotOut

# status message
echo "Beginning analysis..."

# TO-DO: change to formatted round seqs file	
# loop over round sequences
while read seqData; do
	# retrieve the seq
	seq=$(echo $seqData | cut -d"," -f4)
	# update the count data
	countData=$seqData
	# status message
	echo "Processing $seq ..."
	# loop over each round sequences file
	for f2 in $outputsPath"/05_combined/"*_combined\.fa; do
	#for f2 in $outputsPath"/06_formatted/"*_formatted\.fa; do
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
done < $fmtSeqs

# status message
echo "Analysis complete!"
