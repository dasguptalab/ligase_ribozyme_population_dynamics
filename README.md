# RNA Selection Amplification

Project for analyzing RNA sequences from in vitro selection-amplification to isolate RNA ligase ribozyme.

## RNA Evolution Project

The code repository for the analysis pipeline can be found [HERE][1].

### Methods

#### Assignment of sequence families
Sequencing reads from each round were pre-processed using custom BASH and R scripts that are available on GitHub (https://github.com/ElizabethBrooks/RNA_selection_amplification.git). Paired reads were merged using FLASh (v1.2.11) (Magoč & Salzberg, 2011). The merged sequence files were combined with the unmerged forward read files, since reverse reads were low in quality. Read sequences were filtered by quality (AVGQUAL:30) and adapter content was removed using Trimmomatic (v0.39) (Bolger et al., 2014). Sequences were filtered to keep only those that contain the predefined constant up- and down-stream stem regions (Figure 1) and with 40 bp in-between. The filtered sequences were trimmed to the 40-nt randomized region by removing constant regions on both 5′ and 3′ ends. The resulting trimmed sequences were filtered to remove sequences that appeared less than 10 times, then de-replicated to keep only the unique sequences. After pre-processing, between 92,366 and 1,063,996 unique sequences remained per round. Unique sequences from rounds 6 to 8 were clustered into sequence families using Clustal Omega (v1.2.4) (Sievers & Higgins, 2018), which uses mBed-like clustering for a given maximum number of sequences per cluster (N = 500). Since we were interested in tracking conservation or sustained variability of all nucleotide positions across all rounds of selection, sequences from all eight rounds of selection with greater than or equal to  90% similarity to the most abundant peak sequences in the final selection round (round 8) were quantified.

#### Progress of selection across eight rounds (2 March 2025)

| Round 1 | Round 2 | Round 3 | Round 4 | Round 5 | Round 6 | Round 7 | Round 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Reaction Time (min)| 120 | 60 |30 | 20 | 30 | 10 | 10 |10 |
| \[Mg<sup>2+</sup>] (mM) | 20 | 20 | 20 | 20 | 20 | 20 | 10 | 5 |
| Total Raw Reads | 1,485,536 | 1,533,916 | 1,649,680 | 1,436,328 | 1,937,410 | 2,336,945 | 1,229,247 | 1,756,169 |
| High Quality Reads | 1,039,660 | 1,067,585 | 1,033,048 | 866,423 | 981,844 | 916,485 | 582,260 | 889,374 |
| Unique Sequences | 1,036,229 | 1,063,996 | 1,029,483 | 863,123 | 966,495 | 500,507 | 92,366 | 108,529 |
| Percent Diversity | 99.67 | 99.66 | 99.65 | 99.62 | 98.44 | 54.61 | 15.86 | 12.20 |
| Sequence Families | 3019 | 3078 | 2953 | 2549 | 2730 | 1638 | 373 | 402 |

#### Supplemental

##### Analysis Steps

| Step | Analysis | Method |
| --- | --- | --- |
|  1 | Merge the paired-end RNA sequencing reads for each round | FLASh |
|  2 | Combine the merged read files with the unmerged forward read files | BASH |
|  3 | Filter reads by quality and remove any detected adapter content | Trimmomatic |
|  4 | Filter reads by structure to keep only those that contain the expected constant up- and down-stream regions with 40 bp in-between | BASH |
|  5 | Clean reads to retain only the 40 bp in-between region | BASH |
|  6 | Filter to keep unique sequences | BASH |
|  7 | Cluster the unique sequences for each round | Clustal Omega |

### Data Analysis Workflow Steps

The steps 1 through 8 comprise the analysis workflow and correspond to scripts in the code repository. Steps 9 through 11 were performed to prepare data for visualization and plotting.

#### Analysis Scripts

1. merge paired reads for each round using FLASh (see supplement of [this PNAS paper][2])<br>
<b>Note:</b> reads were not merged in the original analysis and only the forward reads were used in the original analysis workflow (see 0010\_qc\_slx.py from the original analysis code)
2. combine merged read files with the unmerged forward read files using BASH<br>
<b>Note:</b> reverse reads are lower quality and typically do not pass filtering by quality or structure.
<b>To-do:</b> update scripts to reflect the workflow step order change.
3. filter reads by quality (AVGQUAL:30) and remove detected adapter content using Trimmomatic<br>
<b>Note:</b> the quality filtering is similar in approach as the original analysis (see 0010\_qc\_slx.py from the original analysis code).
4. filter reads by structure to keep only those that contain the expected constant up- and down-stream regions (with 40 bp in-between in this workflow, but not the original) using BASH ({N}CGGTAGGTCCCTTAGCCAAAAAAGGACAGCG{40}CGCTGTCCGT{M})<br>
<b>Note:</b> the defined start sequence was "ACGGACAGCG" and the end sequence was "CGCTGTCCTTTTTTGGCTAAGGGACCTACCG". Additionally, it is recommended to simply reverse complement the constant regions instead of all reads, at this stage in the analysis (see 0010\_qc\_slx.py from the original analysis code).
5. clean reads to retain only the 40 bp in-between region using BASH<br>
<b>Note:</b> there does not appear to be a filter in the original analysis to keep reads that are only 40-bp in the original analysis workflow.
6. remove sequences that appear less than 10 times (see 0015\_g10\_seqs.py from the original analysis code) and re-format reads and headers using BASH<br>
<b>Note:</b> only unique reeds were kept in this analysis and the sequence headers were updated to contain the run name, arbitrary sequence ID, and read counts for the unique sequence.
7. cluster read sequences for each round (run) using Clustal Omega<br>
	<b>07a.</b>  cluster with a soft maximum of 500 sequences in sub-clusters (cluster-size=500), which is what the original analysis used (see 0020\_cluster\_slxn.py from the original analysis code)<br>
	<b>07b.</b>  cluster with the default soft maximum of 100 sequences in sub-clusters (see [Clustal Omega README][3])
8. create tables with the reverse compliment of cluster sequences and peak sequences, in addition to the cluster and sequence information (run name, sequence ID, read counts, cluster ID, sequence counts, reverse complimented sequence)
9. create tables with the statistics (average, standard deviation, highest, lowest) for the percent identity of cluster sequences relative to the peak sequence within each cluster (see the JAX's [Introduction to Sequence Comparison][4])
10. create tables with counts of the number of sequences in each round 8 cluster sequence family across sequencing rounds
11. create tables with counts of the number of sequences shared across sequencing rounds for:<br>
	<b>11a.</b> all sequences that appear at least 10 times per round<br>
	<b>11b.</b> the top 10 sequences per round

<b>Note</b> that the 00a\_qc.sh script can be used to assess the quality of the fastq data after the 01\_merged, 02\_trimmed, and 03\_filtered stages. The 00b\_analyze.sh script can be used to assess the resulting fasta data from the 01\_merged, 02\_trimmed, 03\_filtered, 04\_cleaned, 05\_combined, and 06\_formatted stages.

### Data Visualization Workflow Steps

The following steps were taken to reproduce tables and plots from the slides (see the "For Elizabeth\_Population dynamics during RNA evolution.pptx" file) and additional interesting tables and plots.

<b>Note:</b> the round 8 cluster sequence families include sequences at least 90 percent identical to the cluster peak (see the 09\_identify.sh script).

#### Plotting Scripts

1. create a line plot with the percent unique sequences (see slide 4)
2. create line plots with the round 8 cluster sequence family counts and abundances (see slide 4)
3. create hetamaps with the round 8 cluster sequence family counts and abundances (expanding on slide 4)
4. create hetamaps with the log counts for the top10 sequences per round (see slide 5)
5. create hetamaps with the log counts for all sequences that appear at least 10 times per round (expanding on slide 5)
6. create hetamaps with the round 8 cluster sequence family base conservation with the substrate 3' overhang bases highlighted (see slide 6)
7. create heatmaps with the base conservation of the substrate 3' overhang for all sequences that appear at least 10 times per round (see slide 7)
8. create line plots with the summed counts and abundances of the top 10 sequences per round (expanding on slide 5)
9. create heatmaps with the round 8 cluster sequence family stem base pairing conservation (see slide 8)
10. create sequence logos with the round 8 cluster family sequences
11. create combined plots with the sequence logos, base conservation heatmap, and stem base pairing conservation heatmaps for the round 8 cluster family sequences

<b>Note:</b> the sequence data analysis plots with the round 8 cluster sequence families were produced from stage 07a. Additionally, the 00a\_cluster\_sequence\_identity.R and 00b\_cluster\_sequence\_identity.R scripts can be used to analyze the percent identities across clusters from stages 07a and 07b respectively.

#### Progress Assessment
For analysis steps 01 to 07 use BASH to:<br>
<b>00a.</b>  create QC reports for each set of raw and processed data using FastQC and MultiQC<br>
<b>00b.</b>  analyze read numbers, read lengths, counts of unique reads, and counts of read names

#### Progress of selection across eight rounds (12 January 2025)

| <div style="width:175px">Statistic</div> | <b>Round 1</b> | <b>Round 2</b> | <b>Round 3</b> | <b>Round 4</b> | <b>Round 5</b> | <b>Round 6</b> | <b>Round 7</b> | <b>Round 8</b> |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Reaction Time (min)| 120 | 60 |30 | 20 | 30 | 10 | 10 |10 |
| \[Mg<sup>2+</sup>] (mM) | 20 | 20 | 20 | 20 | 20 | 20 | 10 | 5 |
| Total Raw Reads | 1,485,536 | 1,533,916 | 1,649,680 | 1,436,328 | 1,937,410 | 2,336,945 | 1,229,247 | 1,756,169 |
| High Quality Reads | 1,039,660 | 1,067,585 | 1,033,048 | 866,423 | 981,844 | 916,485 | 582,260 | 889,374 |
| Unique Reads | 1,036,229 | 1,063,996 | 1,029,483 | 863,123 | 966,495 | 500,507 | 92,366 | 108,529 |
| Diversity (%) | 99.67 | 99.66 | 99.65 | 99.62 | 98.44 | 54.61 | 15.86 | 12.20 |
| Sequences with > 9 Reads | 5 | 3 | 5 | 4 | 283 | 4,001 | 2,703 | 2,100 |
| Sequence Families | NA | NA | NA | NA | NA | 16 | 14 | 13 |

#### Analysis Workflow Progress - Sequence Statistics (28 October 2024)

| | Statistic | Round 1 | Round 2 | Round 3 | Round 4 | Round 5 | Round 6 | Round 7 | Round 8 | Doped 1 | Doped 2 | Doped 3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | Reaction Time (min)| 120 | 60 |30 | 20 | 30 | 10 | 10 |10 | --- | --- | --- |
| B | [Mg<sup>2+</sup>] (mM) | 20 | 20 | 20 | 20 | 20 | 20 | 10 | 5 | --- | --- | --- |
| C | Total Raw Reads | 1,485,536 | 1,533,916 | 1,649,680 | 1,436,328 | 1,937,410 | 2,336,945 | 1,229,247 | 1,756,169 | 1,226,539 | 1,090,909 | 1,656,088 |
| D | High Quality Reads | 1,039,660 | 1,067,585 | 1,033,048 | 866,423 | 981,844 | 916,485 | 582,260 | 889,374 | 865,509 | 807,849 | 1,143,871 |
| E | Unique Reads | 1,036,229 | 1,063,996 | 1,029,483 | 863,123 | 966,495 | 500,507 | 92,366 | 108,529 | 842,149 | 746,445 | 988,626 |
| F | Diversity (%) | 99.67 | 99.66 | 99.65 | 99.62 | 98.44 | 54.61 | 15.86 | 12.20 | 97.30 | 92.40 | 86.43 |
| G | Sequences with > 9 Reads | 5 | 3 | 5 | 4 | 283 | 4,001 | 2,703 | 2,100 | 62 | 749 | 1,690 |
| H | Sequence Families | NA | NA | NA | NA | NA | 16 | 14 | 13 | NA | 4 | 6 |

##### Additional Analysis Results
| | Statistic | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| G | Sequence Families (100) | NA | NA | NA | NA | 7 | 68 | 45 | 41 | NA | 12 | 27 |

#### qc\_statistics\_table - from original code repository data directory

![qc_statistics_table - from original code repository data directory](images/qc_statistics_table.png)

<b>Note:</b> the number of "Sequences with >10 reads" is actually the number of sequences with greater than 9 reads (see 0015\_g10\_seqs.py from the original analysis code).

#### Progress of in vitro selection - from original slides

![Progress of in vitro selection - from original slides](images/Progress_of_in_vitro_selection.png)

##### Table Data Description
- Total: the number of raw read sequences in either (not both) the forward or reverse read file
- Quality: the number of sequences that passed quality filtering and contain the expected start and end patterns 
- Unique: the number of unique sequences that passed quality filtering and contain the expected start and end patterns
- Diversity%: the percentage of "Unique" sequences contained in the "Quality" sequences
- Families(500): the number of clustered sequence families with the original soft maximum of 500 sequences in sub-clusters
- Families(100): the number of clustered sequence families with the default soft maximum of 100 sequences in sub-clusters

### TO-DO
- Add software versions
- Add the reverse complimenting of sequences to step 6 for downstream analysis


[1]: https://github.com/ElizabethBrooks/RNA_selection_amplification
[2]: https://doi.org/10.1073/pnas.2321592121
[3]: https://github.com/GSLBiotech/clustal-omega/blob/master/README
[4]: https://www.jax.org/-/media/jaxweb/files/education-and-learning/ttgg-seq-comparison/seqcomp_introduction.pdf?rev=93c9ff2010234a4bb0dfa0ed043de28e#:~:text=%F0%9D%91%83%F0%9D%91%92%F0%9D%91%9F%F0%9D%91%90%F0%9D%91%92%F0%9D%91%9B%F0%9D%91%A1%20%F0%9D%90%BC%F0%9D%91%91%F0%9D%91%92%F0%9D%91%9B%F0%9D%91%A1%F0%9D%91%96%F0%9D%91%A1%F0%9D%91%A6%20%3D%20%23%20%F0%9D%91%9D%F0%9D%91%9C%F0%9D%91%A0%F0%9D%91%96%F0%9D%91%A1%F0%9D%91%96%F0%9D%91%9C%F0%9D%91%9B%F0%9D%91%A0%20%E2%88%92%20%23,of%20sequences%20can%20be%20compared.
