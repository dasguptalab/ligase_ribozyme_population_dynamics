#!/bin/bash

# script to perform fastqc quality control of paired end reads
# usage: bash 00b_report.sh analysisType
# usage ex: bash 00b_report.sh raw
# usage ex: bash 00b_report.sh trimmed
# usage ex: bash 00b_report.sh merged
# usage ex: bash 00b_report.sh combined

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
	multiqc $qcOut -o $outputsPath
elif [[ $analysisType == "trimmed" ]]; then
	# run multiqc on all
	multiqc $qcOut -o $outputsPath -n "qc_trimmed_all"
	# run multiqc on paired data
	multiqc $qcOut"/"*_p* -o $outputsPath -n "qc_trimmed_paired"
else
	# run multiqc
	multiqc $qcOut -o $outputsPath -n "qc_"$analysisType
fi

#Print status message
echo "Analysis complete!"
