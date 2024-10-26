#!/bin/bash

# script to perform fastqc quality control of paired end reads
# usage: bash 00b_report.sh analysisType
# usage ex: bash 00b_report.sh raw
# usage ex: bash 00b_report.sh merged
# usage ex: bash 00b_report.sh trimmed
# usage ex: bash 00b_report.sh filtered

# retrieve input analysis type
analysisType=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# retrieve directory for analysis
outputsPath=$outputsPath"/qc_"$analysisType

# move to the new directory
cd $outputsPath

# status message
echo "Processing..."

# check input analysis type
if [[ $analysisType == "raw" ]]; then
	# run multiqc
	multiqc $outputsPath -o $outputsPath -n "qc_raw"
elif [[ $analysisType == "trimmed"* ]]; then
	# run multiqc on all
	multiqc $outputsPath -o $outputsPath -n "qc_"$analysisType"_all"
	# run multiqc on paired data
	#multiqc $outputsPath"/"*_p* -o $outputsPath -n "qc_"$analysisType"_paired"
else
	# run multiqc
	multiqc $outputsPath -o $outputsPath -n "qc_"$analysisType
fi

# clean up
#rm -r $outputsPath

#Print status message
echo "Analysis complete!"
