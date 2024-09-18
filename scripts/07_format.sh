#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_format_jobOutput
#$ -q largemem

# script to subset sequences and format headers
# usage: qsub 07_format.sh inputFile
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs_s4q15/filtered_combined/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 07_format.sh "${fileList[$i]}"; done
## job 816500
## job 816501
## job 816502
## job 816503
## job 816504
## job 816505
## job 816506
## job 816507
## job 816508
## job 816509
## job 816511
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/cleaned_s4q20/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 07_format.sh "${fileList[$i]}"; done
## job 819632
## job 819633
## job 819634
## job 819635
## job 819636
## job 819637
## job 819638
## job 819639
## job 819640
## job 819641
## job 819642

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
formatOut=$outputsPath"/formatted_subset"
mkdir $formatOut

# move to the new directory
cd $formatOut

# create sequence run name
#newName=$(basename $f1 | sed "s/_L001_p.*\.flt\.fq//g")
# clean up file name
fileName=$(basename $inputFile | sed "s/\.fa/\.fmt\.fa/g")

# pre-clean up
rm $formatOut"/"$fileName

# status message
echo "Beginning analysis..."

# subset read sequences and re-format headers
while read -r line; do
	# check if the current line is a header
	if [[ $line == "@"* ]]; then # read header
		# save previous line, which should be the leader of the current sequence
		headerLine=$line
	else # read sequence
		# status message
		echo "Processing $line ..."
		# get read count
		readCount=$(cat $inputFile | grep $line | wc -l | tr -d " ")
		# check if there are more than 10 of the current read
		if (( $readCount > 10 )); then # add read info to outputs
			# re-format read headers and add to outputs
			echo $headerLine | cut -d" " -f1 | sed "s/^@/>/g" >> $formatOut"/"$fileName
			# add sequence to outputs
			echo $line >> $formatOut"/"$fileName
		fi
	fi
done < $inputFile

# status message
echo "Analysis complete!"
