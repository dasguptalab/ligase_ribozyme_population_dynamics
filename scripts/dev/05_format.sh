#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_format_jobOutput
#$ -q largemem

# script to subset sequences and format headers
# usage: qsub 05_format.sh inputFile
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs_s4q15/filtered_combined/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 05_format.sh "${fileList[$i]}"; done
## jobs 816500 to 816509 and 816511
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/cleaned_s4q20/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 05_format.sh "${fileList[$i]}"; done
## jobs 819632 to 819642
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/cleaned_merged/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 05_format.sh "${fileList[$i]}"; done
## jobs 870685 to 870695
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/cleaned_merged_copy/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 05_format.sh "${fileList[$i]}"; done
## jobs 871478 to 871490

# retrieve input file
inputFile=$1

# retrieve input analysis type
analysisType=$(basename $1)

# retrieve the analysis type
analysisTag=$(grep "analysis:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/analysis://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# make a new directory for analysis
formatOut=$outputsPath"/formatted_"$analysisTag
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

# very innificient... consider using sort to prepare the data and reduce search time
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
