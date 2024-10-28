# RNA Selection Amplification

Repository for analyzing RNA sequences from in vitro selection-amplification to isolate RNA ligase ribozyme.

## RNA Evolution Project

The code for the analysis pipeline can be found [HERE](https://github.com/ElizabethBrooks/RNA_selection_amplification).

### Analysis Workflow

0. after each analysis step and for each run use BASH to:
- create QC reports for each set of raw and processed data using FastQC and MultiQC
- analyze read numbers and lengths 
1. merge the trimmed paired reads for each run using FLASh (see supplement of https://doi.org/10.1073/pnas.2321592121)
2. filter reads by quality (AVGQUAL:30) and remove detected adapter content using Trimmomatic
3. filter reads by structure to keep only those that contain the expected constant up- and down-stream regions with 40 bp in-between using BASH
4. clean reads to retain only the 40 bp in-between region using BASH
5. combine merged read files with the unmerged forward read files using BASH
6. remove sequences that appear less than 10 times and re-format read headers using BASH
7. cluster read sequences for each run using Clustal Omega
8. reproduce tables and plots from slides/paper using BASH and R
9. create additional tables or plats, as needed

### Progress of in vitro selection (28 October 2024)

| Stat | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Total | 1,485,536 | 1,533,916 | 1,649,680 | 1,436,328 | 1,937,410 | 2,336,945 | 1,229,247 | 1,756,169 | 1,226,539 | 1,090,909 | 1,656,088 |
| Quality | 1,039,660 | 1,067,585 | 1,033,048 | 866,423 | 981,844 | 916,485 | 582,260 | 889,374 | 865,509 | 807,849 | 1,143,871 |
| Unique | 1,036,229 | 1,063,996 | 1,029,483 | 863,123 | 966,495 | 500,507 | 92,366 | 108,529 | 842,149 | 746,445 | 988,626 |
| Diversity% | 99.67 | 99.66 | 99.65 | 99.62 | 98.44 | 54.61 | 15.86 | 12.20 | 97.30 | 92.40 | 86.43 |
| Families | NA | NA | NA | NA | NA | 15 | 13 | 12 | NA | 3 | 5 |

![Progress of in vitro selection - OG](images/Progress_of_in_vitro_selection.png)
