# RNA Selection Amplification

Repository for analyzing RNA sequences from in vitro selection-amplification to isolate RNA ligase ribozyme.

## RNA Evolution Project

The code for the analysis pipeline can be found [HERE](https://github.com/ElizabethBrooks/RNA_selection_amplification).

### Analysis Workflow

1. merge paired reads for each run using FLASh (see supplement of [this PNAS paper](https://doi.org/10.1073/pnas.2321592121))<br>
<b>Note:</b> reads were not merged and only the forward reads were used in the original analysis workflow (see 0010\_qc\_slx.py from the original analysis code)
2. filter reads by quality (AVGQUAL:30) and remove detected adapter content using Trimmomatic<br>
<b>Note:</b> the quality filtering is similar in approach as the original analysis (see 0010\_qc\_slx.py from the original analysis code)
3. filter reads by structure to keep only those that contain the expected constant up- and down-stream regions (with 40 bp in-between) using BASH<br>
<b>Note:</b> why did they take the reverse compliment of all reads rather than simply the constant regions? (see 0010\_qc\_slx.py from the original analysis code)
4. clean reads to retain only the 40 bp in-between region using BASH<br>
<b>Note:</b> I did not see a filter to keep reads that are only 40-bp in the original analysis workflow
5. combine merged read files with the unmerged forward read files using BASH<br>
<b>Note:</b> reverse reads are lower quality and typically do not pass filtering by quality or structure 
6. remove sequences that appear less than 10 times (see 0015\_g10\_seqs.py from the original analysis code) and re-format reads and headers using BASH<br>
<b>Note:</b> only unique reeds were kept and the sequence headers were updated to contain the run number, arbitrary sequence number, and read counts for the unique sequence
7. cluster read sequences for each run using Clustal Omega 
	- cluster with a soft maximum of 500 sequences in sub-clusters (cluster-size=500), which is what the original analysis used (see 0020\_cluster\_slxn.py from the original analysis code)
	- cluster with the default soft maximum of 100 sequences in sub-clusters (see [Clustal Omega README](https://github.com/GSLBiotech/clustal-omega/blob/master/README))
8. for each analysis step and for each run use BASH to:
	- create QC reports for each set of raw and processed data using FastQC and MultiQC
	- analyze read numbers and lengths 
9. reproduce tables and plots from slides/paper using BASH and R (see 0030\_mk\_qc\_table.py from the original analysis code)
10. create additional tables and plots, as needed

### Progress of in vitro selection (28 October 2024)

| Stat | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Total | 1,485,536 | 1,533,916 | 1,649,680 | 1,436,328 | 1,937,410 | 2,336,945 | 1,229,247 | 1,756,169 | 1,226,539 | 1,090,909 | 1,656,088 |
| Quality | 1,039,660 | 1,067,585 | 1,033,048 | 866,423 | 981,844 | 916,485 | 582,260 | 889,374 | 865,509 | 807,849 | 1,143,871 |
| Unique | 1,036,229 | 1,063,996 | 1,029,483 | 863,123 | 966,495 | 500,507 | 92,366 | 108,529 | 842,149 | 746,445 | 988,626 |
| Diversity% | 99.67 | 99.66 | 99.65 | 99.62 | 98.44 | 54.61 | 15.86 | 12.20 | 97.30 | 92.40 | 86.43 |
| Families(500) | NA | NA | NA | NA | NA | 15 | 13 | 12 | NA | 3 | 5 |
| Families(100) | NA | NA | NA | NA | 6 | 67 | 44 | 40 | NA | 11 | 26 |

![Progress of in vitro selection - OG](images/Progress_of_in_vitro_selection.png)

#### Table Data Description
- Total: the number of raw read sequences in either (not both) the forward or reverse read file
- Quality: the number of sequences that passed quality filtering and contain the expected start and end patterns 
- Unique: the number of unique sequences that passed quality filtering and contain the expected start and end patterns
- Diversity%: the percentage of "Unique" sequences contained in the "Quality" sequences
- Families(500): the number of clustered sequence families with the original soft maximum of 500 sequences in sub-clusters
- Families(100): the number of clustered sequence families with the default soft maximum of 100 sequences in sub-clusters
