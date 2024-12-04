# RNA Selection Amplification

Project for analyzing RNA sequences from in vitro selection-amplification to isolate RNA ligase ribozyme.

## RNA Evolution Project

The code repository for the analysis pipeline can be found [HERE][1].

### Data Analysis Workflow

The following steps comprise the analysis workflow and correspond to scripts in the code repository.

#### Analysis

1. merge paired reads for each run using FLASh (see supplement of [this PNAS paper][2])<br>
<b>Note:</b> reads were not merged in the original analysis and only the forward reads were used in the original analysis workflow (see 0010\_qc\_slx.py from the original analysis code)
2. filter reads by quality (AVGQUAL:30) and remove detected adapter content using Trimmomatic<br>
<b>Note:</b> the quality filtering is similar in approach as the original analysis (see 0010\_qc\_slx.py from the original analysis code)
3. filter reads by structure to keep only those that contain the expected constant up- and down-stream regions (with 40 bp in-between in this workflow, but not the original) using BASH<br>
<b>Note:</b> it is recommended to simply reverse complement the constant regions instead of all reads, at this stage in the analysis (see 0010\_qc\_slx.py from the original analysis code)
4. clean reads to retain only the 40 bp in-between region using BASH<br>
<b>Note:</b> there does not appear to be a filter in the original analysis to keep reads that are only 40-bp in the original analysis workflow
5. combine merged read files with the unmerged forward read files using BASH<br>
<b>Note:</b> reverse reads are lower quality and typically do not pass filtering by quality or structure 
6. remove sequences that appear less than 10 times (see 0015\_g10\_seqs.py from the original analysis code) and re-format reads and headers using BASH<br>
<b>Note:</b> only unique reeds were kept in this analysis and the sequence headers were updated to contain the run name, arbitrary sequence ID, and read counts for the unique sequence
7. cluster read sequences for each run using Clustal Omega<br>
	<b>07a.</b>  cluster with a soft maximum of 500 sequences in sub-clusters (cluster-size=500), which is what the original analysis used (see 0020\_cluster\_slxn.py from the original analysis code)<br>
	<b>07b.</b>  cluster with the default soft maximum of 100 sequences in sub-clusters (see [Clustal Omega README][3])
8. create tables with the reverse compliment of cluster sequences and peak sequences, in addition to the cluster and sequence information (run name, sequence ID, read counts, cluster ID, sequence counts, reverse complimented sequence)
9. create tables with the statistics (average, standard deviation, highest, lowest) for the percent identity of cluster sequences relative to the peak sequence within each cluster (see the JAX's [Introduction to Sequence Comparison][4])

<b>Note</b> that the 00a\_qc.sh script can be used to assess the quality of the fastq data after the 01\_merged, 02\_trimmed, and 03\_filtered stages. The 00b\_analyze.sh script can be used to assess the resulting fasta data from the 01\_merged, 02\_trimmed, 03\_filtered, 04\_cleaned, 05\_combined, and 06\_formatted stages.

### Data Visualization Workflow

The following steps reproduce tables and plots from the slides/paper (see 0030\_mk\_qc\_table.py from the original analysis code) and additional interesting tables and plots.

#### Analysis

10. count the number of sequences in each round 8 cluster sequence family across sequencing rounds
11. count the number of sequences shared across sequencing rounds for:<br>
	<b>11a.</b> all sequences that appear at least 10 times per round<br>
	<b>11b.</b> the top 10 sequences per round

<b>Note:</b> the round 8 cluster sequence families include sequences at least 90 percent identical to the cluster peak (see the 09\_identify.sh script).

#### Plotting

1. create a line plot with the percent unique sequences
2. create line plots with the round 8 cluster sequence family fraction abundances and log fraction abundances
3. create hetamaps with the round 8 cluster sequence family log counts and log fraction abundances
4. create hetamaps with the log counts for the top10 sequences per round
5. create hetamaps with the log counts for all sequences that appear at least 10 times per round
6. create hetamaps with the round 8 cluster sequence family base conservation
7. create hetamaps with the round 8 cluster sequence family base conservation of the substrate 3' overhang
8. create hetamaps with the base conservation of the substrate 3' overhang for all sequences that appear at least 10 times per round

<b>Note:</b> the sequence data analysis plots with the round 8 cluster sequence families were produced from stage 07a. Additionally, the 00a\_cluster\_sequence\_identity.R and 00b\_cluster\_sequence\_identity.R scripts can be used to analyze the percent identities across clusters from stages 07a and 07b respectively.

#### Progress Assessment
For analysis steps 01 to 07 use BASH to:<br>
<b>00a.</b>  create QC reports for each set of raw and processed data using FastQC and MultiQC<br>
<b>00b.</b>  analyze read numbers, read lengths, counts of unique reads, and counts of read names

#### Analysis Workflow Progress - Sequence Statistics (28 October 2024)

| Statistic | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Total | 1,485,536 | 1,533,916 | 1,649,680 | 1,436,328 | 1,937,410 | 2,336,945 | 1,229,247 | 1,756,169 | 1,226,539 | 1,090,909 | 1,656,088 |
| Quality | 1,039,660 | 1,067,585 | 1,033,048 | 866,423 | 981,844 | 916,485 | 582,260 | 889,374 | 865,509 | 807,849 | 1,143,871 |
| Unique | 1,036,229 | 1,063,996 | 1,029,483 | 863,123 | 966,495 | 500,507 | 92,366 | 108,529 | 842,149 | 746,445 | 988,626 |
| Diversity% | 99.67 | 99.66 | 99.65 | 99.62 | 98.44 | 54.61 | 15.86 | 12.20 | 97.30 | 92.40 | 86.43 |
| >9Reads | 5 | 3 | 5 | 4 | 283 | 4,001 | 2,703 | 2,100 | 62 | 749 | 1,690 |
| Families(500) | NA | NA | NA | NA | NA | 15 | 13 | 12 | NA | 3 | 5 |
| Families(100) | NA | NA | NA | NA | 6 | 67 | 44 | 40 | NA | 11 | 26 |

#### qc\_statistics\_table - from original code repository data directory

![qc_statistics_table - from original code repository data directory](images/qc_statistics_table.png)

<b>Note:</b> the number of "Sequences with >10 reads" is actually the number of sequences with greater than 9 reads (see 0015\_g10\_seqs.py from the original analysis code)

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
- Create methods summary and references list


[1]: https://github.com/ElizabethBrooks/RNA_selection_amplification
[2]: https://doi.org/10.1073/pnas.2321592121
[3]: https://github.com/GSLBiotech/clustal-omega/blob/master/README
[4]: https://www.jax.org/-/media/jaxweb/files/education-and-learning/ttgg-seq-comparison/seqcomp_introduction.pdf?rev=93c9ff2010234a4bb0dfa0ed043de28e#:~:text=%F0%9D%91%83%F0%9D%91%92%F0%9D%91%9F%F0%9D%91%90%F0%9D%91%92%F0%9D%91%9B%F0%9D%91%A1%20%F0%9D%90%BC%F0%9D%91%91%F0%9D%91%92%F0%9D%91%9B%F0%9D%91%A1%F0%9D%91%96%F0%9D%91%A1%F0%9D%91%A6%20%3D%20%23%20%F0%9D%91%9D%F0%9D%91%9C%F0%9D%91%A0%F0%9D%91%96%F0%9D%91%A1%F0%9D%91%96%F0%9D%91%9C%F0%9D%91%9B%F0%9D%91%A0%20%E2%88%92%20%23,of%20sequences%20can%20be%20compared.
