#!/bin/bash
#$ -r n
#$ -N RNA_randomize_jobOutput
#$ -q largemem

# script to create files with randomized nucelotide sequences
# usage: qsub 14_randomize.sh numNum
# usage example: qsub 14_randomize.sh 1
## jobs 1728185 to 1728188

# retrieve input run num
runNum=$1

# retrieve analysis outputs absolute path
#outputsPath="/Users/bamflappy/PfrenderLab/RNA_evolution/outputs"
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/14_randomized_sequences"

# create outputs directory
mkdir $outDir

# set of nucleotides 
nucSet="ACTG"

# length of seqeuences
seqLength=40

# length of file
numSeqs=1500000

# name out output file
outFile=$outDir"/random_run"$runNum"_sequences.fa"

# pre-clean up
rm $outFile

# loop over each line of the file
for seqNum in $(seq 1 $numSeqs); do
	# status message
	echo "Processing sequence number $seqNum"
	# reset the current sequence
	currSeq=""
	# loop over each nucleotide
	for nucNum in $(seq 1 $seqLength); do
		# generate a random number
		nucIndex=$((RANDOM % 4))
		# select a random nucleotide
		randomNuc="${nucSet:$nucIndex:1}"
		# add the random nucelotide to the sequence
		currSeq="$currSeq$randomNuc"
	done
	# add the randomized sequence to the file
	echo "@random_$seqNum" >> $outFile
	echo $currSeq >> $outFile
done

# status message
echo "Analysis complete!"
