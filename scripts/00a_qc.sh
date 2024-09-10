#!/bin/bash

# script to perform fastqc quality control of paired end reads
# usage: bash 00a_qc.sh analysisType
# usage ex: bash 00a_qc.sh raw
# usage ex: bash 00a_qc.sh trimmed
# usage ex: bash 00a_qc.sh filtered

# retrieve input analysis type
analysisType=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# check input analysis type
if [[ $analysisType == "raw" ]]; then
	# retrieve raw paired reads absolute path for alignment
	readPath=$(grep "pairedReads:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/pairedReads://g")
elif [[ $analysisType == "trimmed" ]]; then
	# retrieve trimmed reads path
	readPath=$outputsPath"/trimmed"
elif [[ $analysisType == "filtered" ]]; then
	# retrieve filtered reads path
	readPath=$outputsPath"/filtered"
fi

# create output results directory
#mkdir $outputsPath

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
echo "Processing..."

# perform QC
fastqc $readPath"/"* -o $qcOut

#Print status message
echo "Analysis complete!"
