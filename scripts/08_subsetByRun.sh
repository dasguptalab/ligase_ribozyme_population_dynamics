#!/bin/bash

# script to filter sequences by length
# usage: bash 08_subsetByRun.sh 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set directory for inputs
formatOut=$outputsPath"/formatted"

# status message
echo "Keeping doped..."

# filter to keep doped reads
# > 36 residues
#cat $formatOut"/combined.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "doped" | sed "s/~/\n/g" > $formatOut"/combined_doped.flt.fmt.fasta"
# 40 residues
#cat $formatOut"/combined.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "doped" | sed "s/~/\n/g" > $formatOut"/combined_doped.flt40.fmt.fasta"

# status message
echo "Removing doped..."

# filter to remove doped reads
# > 36 residues
#cat $formatOut"/combined.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep -v "doped" | sed "s/~/\n/g" > $formatOut"/combined_noDoped.flt.fmt.fasta"
# 40 residues
#cat $formatOut"/combined.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep -v "doped" | sed "s/~/\n/g" > $formatOut"/combined_noDoped.flt40.fmt.fasta"

# status message
echo "Subsetting by run..."

# filter to subset by run
# > 36 residues
# doped r1
cat $formatOut"/combined_doped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r1_" | sed "s/~/\n/g" > $formatOut"/combined_doped_r1.flt.fmt.fasta"
# doped r2
cat $formatOut"/combined_doped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r2_" | sed "s/~/\n/g" > $formatOut"/combined_doped_r2.flt.fmt.fasta"
# doped r3
cat $formatOut"/combined_doped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r3_" | sed "s/~/\n/g" > $formatOut"/combined_doped_r3.flt.fmt.fasta"
# no doped r1
cat $formatOut"/combined_noDoped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r1_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r1.flt.fmt.fasta"
# no doped r2
cat $formatOut"/combined_noDoped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r2_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r2.flt.fmt.fasta"
# no doped r3
cat $formatOut"/combined_noDoped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r3_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r3.flt.fmt.fasta"
# no doped r4
cat $formatOut"/combined_noDoped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r4_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r4.flt.fmt.fasta"
# no doped r5
cat $formatOut"/combined_noDoped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r5_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r5.flt.fmt.fasta"
# no doped r6
cat $formatOut"/combined_noDoped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r6_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r6.flt.fmt.fasta"
# no doped r7
cat $formatOut"/combined_noDoped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r7_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r7.flt.fmt.fasta"
# no doped r8
cat $formatOut"/combined_noDoped.flt.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r8_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r8.flt.fmt.fasta"
# 40 residues
# doped r1
cat $formatOut"/combined_doped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r1_" | sed "s/~/\n/g" > $formatOut"/combined_doped_r1.flt40.fmt.fasta"
# doped r2
cat $formatOut"/combined_doped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r2_" | sed "s/~/\n/g" > $formatOut"/combined_doped_r2.flt40.fmt.fasta"
# doped r3
cat $formatOut"/combined_doped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r3_" | sed "s/~/\n/g" > $formatOut"/combined_doped_r3.flt40.fmt.fasta"
# no doped r1
cat $formatOut"/combined_noDoped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r1_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r1.flt40.fmt.fasta"
# no doped r2
cat $formatOut"/combined_noDoped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r2_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r2.flt40.fmt.fasta"
# no doped r3
cat $formatOut"/combined_noDoped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r3_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r3.flt40.fmt.fasta"
# no doped r4
cat $formatOut"/combined_noDoped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r4_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r4.flt40.fmt.fasta"
# no doped r5
cat $formatOut"/combined_noDoped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r5_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r5.flt40.fmt.fasta"
# no doped r6
cat $formatOut"/combined_noDoped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r6_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r6.flt40.fmt.fasta"
# no doped r7
cat $formatOut"/combined_noDoped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r7_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r7.flt40.fmt.fasta"
# no doped r8
cat $formatOut"/combined_noDoped.flt40.fmt.fasta" | tr '\n' '~' | sed "s/~>/\n>/g" | grep "r8_" | sed "s/~/\n/g" > $formatOut"/combined_noDoped_r8.flt40.fmt.fasta"

# status message
echo "Analysis complete!"
