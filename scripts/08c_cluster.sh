#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N RNA_cluster_maxnumseq_clustalo_jobOutput
#$ -pe smp 8
#$ -q largemem

# script to cluster sequences using clustalo and --maxnumseq=1300000
# usage: qsub 08c_cluster.sh inputFile
# usage ex: for i in /scratch365/ebrooks5/RNA_evolution/outputs/formatted/*; do qsub 08c_cluster.sh $i; done
# usage ex: fileList=(/scratch365/ebrooks5/RNA_evolution/outputs/formatted/*); for ((i=${#fileList[@]}-1; i>=0; i--)); do qsub 08c_cluster.sh "${fileList[$i]}"; done
## job 814075 -> FATAL: Memory allocation for distance matrix failed
## job 814076 -> FATAL: Memory allocation for distance matrix failed
## job 814077 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_maxnumseq_1300000_r6_S6_L001_combined_fmt/clustered_r6_S6_L001_combined_fmt
## job 814078
## job 814079 -> /opt/sge/crc/spool/d32cepyc245/job_scripts/814079: line 41: 1293623 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --maxnumseq=1300000
## job 814081 -> /opt/sge/crc/spool/d32cepyc229/job_scripts/814081: line 41: 2429763 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --maxnumseq=1300000
## job 814082 -> /opt/sge/crc/spool/d32cepyc207/job_scripts/814082: line 41: 2186459 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --maxnumseq=1300000
## job 814083 -> /opt/sge/crc/spool/d32cepyc193/job_scripts/814083: line 41: 2886212 Killed                  clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --maxnumseq=1300000
## job 814084 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_maxnumseq_1300000_doped21-r3_S12_L001_combined_fmt/clustered_doped21-r3_S12_L001_combined_fmt
## job 814085 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_maxnumseq_1300000_doped21-r2_S11_L001_combined_fmt/clustered_doped21-r2_S11_L001_combined_fmt
## job 814086 -> Alignment written to /scratch365/ebrooks5/RNA_evolution/outputs/clustered_maxnumseq_1300000_doped21-r1_S10_L001_combined_fmt/clustered_doped21-r1_S10_L001_combined_fmt

# load the software module
module load bio/0724

# retrieve input file
inputFile=$1

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# clean up input file name
nameTag=$(basename $inputFile | sed "s/\.fa//g" | sed "s/\./_/g")

# retrieve the analysis type
analysisType=$(dirname $inputFile)
analysisType=$(basename $clusterOut)

# make a directory for the clustering outputs
clusterOut=$outputsPath"/clustered_maxnumseq_1300000_"$analysisType
mkdir $clusterOut

# make a new directory for analysis
clusterOut=$clusterOut"/"$nameTag
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
clustalo --threads=$NSLOTS -v -i $inputFile -o $clusterOut"/clustered_"$nameTag --maxnumseq=1300000

# status message
echo "Analysis complete!"
