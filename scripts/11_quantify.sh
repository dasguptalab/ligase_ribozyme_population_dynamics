#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_quantify_jobOutput
#$ -q largemem

# script to filter fastq files and keep sequences with matching up- and down-stream sequences
# usage: qsub 11_quantify.sh inputFile inputRun
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_formatted.fa; do runInput=$(basename $i | sed "s/_formatted.fa//g"); qsub 11_quantify.sh 07a_clustered $runInput; done
## job 1004471 to 1004481
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_formatted.fa; do runInput=$(basename $i | sed "s/_formatted.fa//g"); qsub 11_quantify.sh 07b_clustered $runInput; done
## job 1004483 to 1004493
# alternate usage ex: for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/*_formatted.fa; do runInput=$(basename $i | sed "s/_formatted.fa//g"); bash 11_quantify.sh 07a_clustered $runInput; done
# alternate usage ex: for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/06_formatted/*_formatted.fa; do runInput=$(basename $i | sed "s/_formatted.fa//g"); bash 11_quantify.sh 07b_clustered $runInput; done

# retrieve input file
inputFile=$1

# retrieve input run name
inputRun=$2

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")
#analysisTag=$(grep "analysis:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve the inputs path
inputsPath=$outputsPath"/06_formatted"

# retrieve input sequences
inputSeqs=$inputsPath"/"$inputRun"_formatted.fa"

# name of a new directory for analysis
tablesOut=$outputsPath"/10_quantified"
# process just the top 10 most abundant sequences
#tablesOut=$outputsPath"/10_quantified_top10"

# make a new directory for analysis
mkdir $tablesOut

# make a new directory for analysis
tablesOut=$tablesOut"/"$inputFile
mkdir $tablesOut

# move to outputs directory
cd $tablesOut

# name formatted sequences file
fmtSeqs=$inputsPath"/"$inputRun"_formatted.tmp.fa"

# re-format input sequences for processing
cat $inputSeqs | tr "\n" "," | sed "s/>/\n>/g" | sed "s/,$//g" | sed '/^[[:space:]]*$/d' > $fmtSeqs
# process just the top 10 most abundant sequences
#cat $inputSeqs | tr "\n" "," | sed "s/>/\n>/g" | sed "s/,$//g" | sed '/^[[:space:]]*$/d' | sort -k3 -n | head -10 > $fmtSeqs

# add final new line
#echo "" >> $fmtSeqs

# name output file
countsOut=$tablesOut"/"$inputRun"_counts_table.csv"
countsPlotOut=$tablesOut"/"$inputRun"_counts_plot_table.csv"

# retrieve header
inputHeader="run_name,sequence_ID,read_counts,sequence"

# add the run names to the header 
header=$(echo $inputHeader",doped21-r1_counts,doped21-r2_counts,doped21-r3_counts,r1_counts,r2_counts,r3_counts,r4_counts,r5_counts,r6_counts,r7_counts,r8_counts")
headerPlot=$(echo $inputHeader",counts,counts_run_name")

# add a header to the counts data tmp file
echo $header > $countsOut
echo $headerPlot > $countsPlotOut

# status message
echo "Beginning analysis of $inputRun ..."

# TO-DO: change to formatted round seqs file	
# loop over round sequences
while read data; do
	# clean up the run name
	seqData=$(echo $data | sed "s/>r//g")
	# remove the sequence
	dataNoSeq=$(echo $seqData | cut -d"," -f1-3)
	# retrieve the seq
	seq=$(echo $seqData | cut -d"," -f4)
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
done < $fmtSeqs

# add run tags to sequence IDs
#for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_quantified_top10/07a_clustered/r*_counts_plot_table.csv; do runTag=$(basename $i | cut -d"_" -f1 | sed "s/r//g"); tail -n+2 $i | awk -v runIn=$runTag 'BEGIN{FS=OFS=","}{$2 = runIn"_"$2; print}' > $i.fmt; done

# after processing the last round of data, combine all plotting data files
#head -1 /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_quantified_top10/07a_clustered/r8_S8_L001_counts_plot_table.csv > /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_quantified_top10/07a_clustered/counts_plot_table_noDoped.csv
#for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_quantified_top10/07a_clustered/r*_counts_plot_table.csv.fmt; do tail -n+2 $i | grep -v "doped" >> /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_quantified_top10/07a_clustered/counts_plot_table_noDoped.csv; done

# clean up
#rm /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10_quantified_top10/07a_clustered/r*_counts_plot_table.csv.fmt

# status message
echo "Analysis complete!"
