#!/bin/bash
#$ -r n
#$ -N RNA_quantify_a_jobOutput
#$ -q largemem

# script to count the number of sequences shared across runs
# usage: qsub 09a_quantify.sh inputRun runName
# usage ex: bash 09a_quantify.sh r8_S8_L001 r8_S8_L001
## quantification of all sequencess
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/06_formatted/*_formatted_above2.fa; do runInput=$(basename $i | sed "s/_formatted_above2\.fa//g"); echo $runInput; qsub 09a_quantify.sh $runInput; done
## jobs 1273266 to 1273292
## job 1273266 -> doped21-r1_S10_L001 -> 10:02:45:26
## job 1273269 -> doped21-r2_S11_L001 -> 9:00:13:10
## job 1273272 -> doped21-r3_S12_L001 -> 11:19:46:49
## job 1273275 -> r1_S1_L001 -> 12:01:53:05
## job 1273278 -> r2_S2_L001 -> 12:09:43:07
## job 1273281 -> r3_S3_L001 -> 12:01:09:05
## job 1273284 -> r4_S4_L001 -> 10:04:20:18
## job 1273287 -> r5_S5_L001 -> 11:07:11:50
## job 1273289 -> r6_S6_L001 -> 6:01:08:18
## job 1273290 -> r7_S7_L001 -> 1:02:34:23
## job 1273292 -> r8_S8_L001 -> 1:07:36:44

# retrieve input run name
inputRun=$1

# retrieve input run name
runName=$2

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")
#analysisTag=$(grep "analysis:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# setup inputs run data
inputData=$outputsPath"/05_combined/"$runName"_combined.RC.fa"

# retrieve input sequences
inputSeqs=$outputsPath"/06_formatted/"$inputRun"_formatted.fa" ## quantification of all sequencess
#inputSeqs=$outputsPath"/06_formatted/"$inputRun"_formatted_above2.fa"

# name of a new directory for analysis
tablesOut=$outputsPath"/09a_quantified_all"

# make a new directory for analysis
mkdir $tablesOut

# setup tmp inputs data
inputRunData=$tablesOut"/"$inputRun"_"$runName"_combined.RC.tmp.fa"
cat $inputData > $inputRunData

# move to outputs directory
cd $tablesOut

# name formatted sequences file
fmtSeqs=$tablesOut"/"$inputRun"_"$runName"_formatted.tmp.fa"

# re-format input sequences for processing
cat $inputSeqs | tr "\n" "," | sed "s/>/\n>/g" | sed "s/,$//g" | sed '/^[[:space:]]*$/d' > $fmtSeqs

# add final new line
echo "" >> $fmtSeqs

# name output file
#countsOut=$tablesOut"/"$inputRun"_in_"$runName"_counts_table.csv"
countsPlotOut=$tablesOut"/"$inputRun"_in_"$runName"_counts_plot_table.csv"

# retrieve header
inputHeader="run_name,sequence_ID,read_counts,sequence"

# add the run names to the header 
#header=$(echo $inputHeader",doped21-r1_counts,doped21-r2_counts,doped21-r3_counts,r1_counts,r2_counts,r3_counts,r4_counts,r5_counts,r6_counts,r7_counts,r8_counts")
headerPlot=$(echo $inputHeader",counts,counts_run_name")

# add a header to the counts data outputs files
#echo $header > $countsOut
echo $headerPlot > $countsPlotOut

# status message
echo "Beginning analysis of $inputSeqs over $inputRunData ..."

# loop over round sequences
while read data; do
	# clean up the run name
	seqData=$(echo $data | sed "s/>r//g")
	# remove the sequence
	dataNoSeq=$(echo $seqData | cut -d"," -f1-3)
	# retrieve the seq
	seq=$(echo $seqData | cut -d"," -f4)
	# reverse compliment the sequence
	#seqRev=$(echo $seq | tr ACGTacgt TGCAtgca | rev)
	# update the count data
	#countData=$dataNoSeq","$seqRev
	countData=$dataNoSeq","$seq
	countDataOut=$countData
	# status message
	echo "Processing $seq ..."
	# count the number of seq occurances in each round
	numReads=$(cat $inputRunData | grep -wc $seq)
	# add the number of seqs for the round
	countDataOut=$(echo $countDataOut","$numReads)
	# add the counts data to the outputs file
	echo $countData","$numReads","$runName >> $countsPlotOut
	# add the counts data to the outputs file
	#echo $countDataOut >> $countsOut
done < $fmtSeqs

# clean up
rm $fmtSeqs
rm $inputRunData

# add run tags to sequence IDs
#for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/r*_counts_plot_table.csv; do runTag=$(basename $i | cut -d"_" -f1 | sed "s/r//g"); tail -n+2 $i | awk -v runIn=$runTag 'BEGIN{FS=OFS=","}{$2 = runIn"_"$2; print}' > $i.fmt; done
# after processing the last round of data, combine all plotting data files
#head -1 /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/r8_S8_L001_counts_plot_table.csv > /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv
#for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/r*_counts_plot_table.csv.fmt; do tail -n+2 $i | grep -v "doped" >> /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv; done
# clean up
#rm /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/r*_counts_plot_table.csv.fmt

# status message
echo "Analysis complete!"
