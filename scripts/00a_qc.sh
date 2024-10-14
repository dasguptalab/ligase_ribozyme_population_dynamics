#!/bin/bash

# script to perform fastqc quality control of paired end reads
# usage: bash 00a_qc.sh analysisType
# usage ex: bash 00a_qc.sh raw
# usage ex: bash 00a_qc.sh trimmed_s4q20
# usage ex: bash 00a_qc.sh combined_s4q20
# usage ex: bash 00a_qc.sh filtered_s4q20
# usage ex: bash 00a_qc.sh trimmed_merged
# usage ex: bash 00a_qc.sh combined_merged
# usage ex: bash 00a_qc.sh filtered_merged

# load software module for remote servers
module load bio/2.0

# retrieve input analysis type
analysisType=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# check input analysis type
if [[ $analysisType == "raw" ]]; then
	# retrieve raw paired reads absolute path for alignment
	readPath=$(grep "pairedReads:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/pairedReads://g")
else 
	# set reads path
	readPath=$outputsPath"/"$analysisType
fi

#Make a new directory for analysis
qcOut=$outputsPath"/qc_"$analysisType
mkdir $qcOut
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $qcOut directory already exsists... please remove before proceeding."
	exit 1
fi

#Move to the new directory
cd $qcOut

# status message
echo "Beginning analysis..."

# perform QC
fastqc $readPath"/"*\.f*q -o $qcOut

#Print status message
echo "Analysis complete!"
