#!/bin/bash

# script to perform fastqc quality control of paired end reads
# usage: bash 02_multiqc.sh analysisType
# usage ex: bash 02_multiqc.sh raw
# usage ex: bash 02_multiqc.sh trimmed
# usage ex: bash 02_multiqc.sh filtered

# retrieve input analysis type
analysisType=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# retrieve directory for analysis
qcOut=$outputsPath"/qc_"$analysisType

# move to the new directory
cd $qcOut

# status message
echo "Processing..."

# check input analysis type
if [[ $analysisType == "raw" ]]; then
	# run multiqc
	multiqc $qcOut -o $qcOut
elif [[ $analysisType == "trimmed" ]]; then
	# run multiqc on all
	multiqc $qcOut -o $qcOut -n "qc_all"
	# run multiqc on paired data
	multiqc $qcOut"/"*_p* -o $qcOut -n "qc_paired"
elif [[ $analysisType == "filtered" ]]; then
	# run multiqc
	multiqc $qcOut -o $qcOut
fi

#Print status message
echo "Analysis complete!"
