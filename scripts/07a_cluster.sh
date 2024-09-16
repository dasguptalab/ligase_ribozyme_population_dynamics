#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_clustalo_jobOutput
#$ -pe smp 63
#$ -q largemem

# script to cluster sequences using clustalo
# usage: qsub 07a_cluster.sh inputFile
# usage ex: qsub 07a_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs/formatted/r8_S8_L001_combined.fmt.fa
## job 816015 -> FATAL: Memory allocation for distance matrix failed
# usage ex: qsub 07a_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs/formatted/r7_S7_L001_combined.fmt.fa
## job 816016 -> FATAL: Memory allocation for distance matrix failed
# usage ex: qsub 07a_cluster.sh /scratch365/ebrooks5/RNA_evolution/outputs/formatted/r6_S6_L001_combined.fmt.fa
## job 816017
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/formatted/*; do qsub 07a_cluster.sh $i; done
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 07a_cluster.sh "${fileList[$i]}"; done
## job 814037 -> Couldn't allocate MPI memory
## job 814038 -> Couldn't allocate MPI memory
## job 814039 -> /opt/sge/crc/spool/d32cepyc207/job_scripts/814039: line 41: 2143111 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag
## job 814040 
## job 814041 -> /opt/sge/crc/spool/d32cepyc189/job_scripts/814041: line 41: 3944610 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag
## job 814042 -> /opt/sge/crc/spool/d32cepyc214/job_scripts/814042: line 41: 2583840 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag
## job 814043 
## job 814044 -> /opt/sge/crc/spool/d32cepyc217/job_scripts/814044: line 41: 81853 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag
## job 814045 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_doped21-r3_S12_L001_combined_fmt/clustered_doped21-r3_S12_L001_combined_fmt
## job 814046 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_doped21-r2_S11_L001_combined_fmt/clustered_doped21-r2_S11_L001_combined_fmt
## job 814047 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_doped21-r1_S10_L001_combined_fmt/clustered_doped21-r1_S10_L001_combined_fmt

# load the software module
module load bio/0724

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# clean up input file name
nameTag=$(basename $inputFile | sed "s/\.fa//g" | sed "s/\./_/g")

# make a new directory for analysis
clusterOut=$outputsPath"/clustered_"$nameTag
mkdir $clusterOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $clusterOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $clusterOut

# status message
echo "Beginning analysis of $nameTag ..."

# filter to keep sequences with matching up- and down-stream sequences
clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag

# status message
echo "Analysis complete!"
