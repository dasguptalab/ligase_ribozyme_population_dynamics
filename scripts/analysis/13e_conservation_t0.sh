#!/bin/bash
#$ -r n
#$ -N RNA_conserved_t0_jobOutput
#$ -q largemem

# script to run job scripts that identify conserved regions
# usage: qsub 13e_conservation_t0.sh runNum
# usage ex: qsub 13e_conservation_t0.sh 1
## jobs 1732347 to 1732447
## test_23May2025
## jobs 1733954 to 1733957

# load the software
module load bio/0724

# retrieve input run number
runNum=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")
#outputsPath=$(grep "outputs:" ../../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set outputs directory
outDir=$outputsPath"/13e_conservation_t0_run"$runNum
#outDir=$outputsPath"/test_23May2025/13e_conservation_t0_run"$runNum

# create outputs directory
mkdir $outDir

# status message
echo "Processing random$runNum sequences..."

# read in sequence count data for the specified round
seqsFile="/scratch365/ebrooks5/RNA_evolution/outputs/14_randomized_sequences/random"$runNum"_sequences_combined.RC.fa"

# run the analysis
Rscript 13_sequence_conservation.R $outDir $seqsFile

# status message
echo "Analysis complete!"
