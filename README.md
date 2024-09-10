# RNA_selection_amplification

Repository for analyzing RNA sequences from in vitro selection-amplification to isolate RNA ligase ribozyme.

## Analysis Workflow

1. create QC reports for the raw data using FastQC and MultiQC
2. remove low quality reads and the first 23 nucleotides using Trimmomatic
3. merge paired end reads

3. filter reads to keep only those that contain the constant up- (GGACAGCG) and down-stream (CGCTGTCC) regions using BASH
4. remove the constant up- and down-stream regions using BASH
5. create QC reports for the filtered data using FastQC and MultiQC
6. filter to keep only a single read per sequence, preferably the forward read
7. TO-DO: filter reads by length to keep sequences that are 40bp long

8. analyze read numbers and lengths for each run
9. reproduce "Progress of in vitro selection"
- total sequences
- high quality sequences
- unique sequences
- diversity (%)
- number of sequence families
10. reproduce "Different classes of ribozymes after clustering" from the run 8 data

## Run Stats

| Stat | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RawForward | 1485536 | 1533916 | 1649680 | 1436328 | 1937410 | 2336945 | 1229247 | 1756169 | 1226539 | 1090909 | 1656088 |
| RawReverse | 1485536 | 1533916 | 1649680 | 1436328 | 1937410 | 2336945 | 1229247 | 1756169 | 1226539 | 1090909 | 1656088 |
| TrimmedForward | 1125487 | 1177286 | 1267375 | 1107470 | 1482964 | 1800051 | 945991 | 1348014 | 913723 | 812516 | 1241307 |
| TrimmedReverse | 1125487 | 1177286 | 1267375 | 1107470 | 1482964 | 1800051 | 945991 | 1348014 | 913723 | 812516 | 1241307 |
| RegionForward | 838868 | 871694 | 839359 | 728738 | 813261 | 757042 | 479992 | 736784 | 721582 | 671629 | 960867 |
| RegionReverse | 852879 | 881806 | 848978 | 743983 | 822940 | 766913 | 486182 | 748219 | 712159 | 658461 | 951936 |
| UniqueForward | 836538 | 869219 | 836906 | 726372 | 799992 | 378395 | 38617 | 28478 | 705348 | 621216 | 830676 |
| UniqueReverse | 850641 | 879552 | 846795 | 741772 | 810513 | 390711 | 46589 | 35600 | 696551 | 611718 | 828204 |
| Combined | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |

## Progress of in vitro selection

### Forward Paird Reads - Filtered by Region Length

| Stat | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Total | 1485536 | 1533916 | 1649680 | 1436328 | 1937410 | 2336945 | 1229247 | 1756169 | 1226539 | 1090909 | 1656088 |
| Quality | 1125487 | 1177286 | 1267375 | 1107470 | 1482964 | 1800051 | 945991 | 1348014 | 913723 | 812516 | 1241307 |
| Unique | 836538 | 869219 | 836906 | 726372 | 799992 | 378395 | 38617 | 28478 | 705348 | 621216 | 830676 |
| Diversity% | 74 | 74 | 66 | 66 | 54 | 21 | 4 | 2 | 77 | 77 | 67 |
| Families | N/A | N/A | N/A | N/A | N/A | R6 | R7 | R8 | N/A | N/A | N/A |
