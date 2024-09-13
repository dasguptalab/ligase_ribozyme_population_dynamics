# RNA_selection_amplification

Repository for analyzing RNA sequences from in vitro selection-amplification to isolate RNA ligase ribozyme.

## Analysis Workflow

0. after each analysis step and for each run:
- create QC reports for each set of raw and processed data using FastQC and MultiQC
- analyze read numbers and lengths 
1. trim reads to filter by quality reads using Trimmomatic
2. merge the trimmed paired reads for each run using NGmerge
3. filter trimmed unpaired forward and reverse, merged paired, failed merged forward and reverse reads by structure to keep only those that contain the constant up- (GGACAGCG) and down-stream (CGCTGTCC) regions with 40 bp in-between, then retain only the 40 bp in-between region using BASH
4. combine filtered read files with trimmed unpaired reads and make sure there are no duplicate reads using BASH
5. reproduce "Progress of in vitro selection"
- total sequences
- high quality sequences
- unique sequences
- diversity (%)
- number of sequence families
6. reproduce "Different classes of ribozymes after clustering" from the run 8 data

## Progress of in vitro selection - Update (12 September 2024)

To count read names:
for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/\*/\*\.fq; do echo \$i; cat $i | awk 'NR%4==1' | wc -l; done
for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/\*/\*\.fa; do echo \$i; cat $i | awk 'NR%2==1' | wc -l; done

To filter sequences by structure and count read names:
for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/combined/\*\.fq; do echo \$i; cat \$i | grep -Ex -B1 '.\*GGACAGCG.{40}CGCTGTCC.\*' | sed "s/^.\*GGACAGCG//g" | sed "s/CGCTGTCC.\*\$//g" | grep -Ex -B1 '.{40}' | grep -v "^--$" | awk 'NR%2==1 | wc -l; done

## Run Stats

| Stat | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RawForward | 1485536 | 1533916 | 1649680 | 1436328 | 1937410 | 2336945 | 1229247 | 1756169 | 1226539 | 1090909 | 1656088 |
| RawReverse | 1485536 | 1533916 | 1649680 | 1436328 | 1937410 | 2336945 | 1229247 | 1756169 | 1226539 | 1090909 | 1656088 |
| TrimmedPF | 1484182 | 1532658 | 1648235 | 1435225 | 1935700 | 2334968 | 1228211 | 1754397 | 1225204 | 1089533 | 1654299 |
| TrimmedPR | 1484182 | 1532658 | 1648235 | 1435225 | 1935700 | 2334968 | 1228211 | 1754397 | 1225204 | 1089533 | 1654299 |
| TrimmedUF | 1285 | 1196 | 1394 | 1050 | 1632 | 1884 | 976 | 1697 | 1297 | 1352 | 1738 |
| FilteredTUF | 113 | 116 | 207 | 140 | 158 | 113 | 86 | 163 | 186 | 346 | 229 |
| TrimmedUR | 60 | 60 | 48 | 51 | 70 | 85 | 54 | 71 | 38 | 21 | 48 |
| FilteredTUR | 11 | 11 | 12 | 7 | 8 | 24 | 6 | 11 | 5 | 2 | 4 |
| Merged | 913620 | 987397 | 1025706 | 951001 | 1204610 | 1472342 | 799622 | 1049155 | 663634 | 619174 | 911598 |
| FilteredM | 794411 | 824063 | 793456 | 693163 | 767847 | 715212 | 453821 | 697759 | 606866 | 561247 | 810341 |
| UnMergedF | 570562 | 545261 | 622529 | 484224 | 731090 | 862626 | 428589 | 705242 | 561570 | 470359 | 742701 |
| FilteredUMF | 122141 | 128749 | 124080 | 107368 | 125311 | 112480 | 75136 | 115697 | 187021 | 179105 | 246300 |
| UnMergedR | 570562 | 545261 | 622529 | 484224 | 731090 | 862626 | 428589 | 705242 | 561570 | 470359 | 742701 |
| FilteredUMR | 119435 | 118315 | 113713 | 101266 | 109282 | 103501 | 63748 | 100765 | 154110 | 141484 | 205815 |
| FilteredCombined | 1036111 | 1071254 | 1031468 | 901944 | 1002606 | 931330 | 592797 | 914395 | 948188 | 882184 | 1262689 |
| Combined | 2056029 | 2079115 | 2272158 | 1920499 | 2668422 | 3199478 | 1657776 | 2461336 | 1788071 | 1561244 | 2398738 |
| CombinedFiltered | 1036100 | 1071243 | 1031456 | 901937 | 1002598 | 931306 | 592791 | 914384 | 948183 | 882182 | 1262685 |

### Forward Paird Reads - Filtered by Structure

To count unique read sequences:
for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/\*/\*\.fq; do echo \$i; cat $i | awk 'NR%4==0' | sort -u | wc -l; done
for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/\*/\*\.fa; do echo \$i; cat $i | awk 'NR%2==0' | sort -u | wc -l; done

To count unique read names:
for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/\*/\*\.fq; do echo \$i; cat $i | awk 'NR%4==1' | cut -d' ' -f1 | sort -u | wc -l; done
for i in /Users/bamflappy/PfrenderLab/RNA_evolution/outputs/\*/\*\.fa; do echo \$i; cat $i | awk 'NR%2==1' | cut -d' ' -f1 | sort -u | wc -l; done

| Stat | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Total | 1,485,536 | 1,533,916 | 1,649,680 | 1,436,328 | 1,937,410 | 2,336,945 | 1,229,247 | 1,756,169 | 1,226,539 | 1,090,909 | 1,656,088 |
| Quality | 1,036,111 | 1,071,254 | 1,031,468 | 901,944 | 1,002,606 | 931,330 | 592,797 | 914,395 | 948,188 | 882,184 | 1,262,689 |
| Unique | 1,033,568 | 1,068,568 | 1,028,855 | 899,424 | 988,146 | 503,256 | 82,177 | 74,546 | 929,757 | 824,867 | 1,114,276 |
| Diversity% | 100 | 100 | 100 | 100 | 99 | 54 | 14 | 0.1 | 98 | 94 | 88 |
| Families | N/A | N/A | N/A | N/A | N/A | R6 | R7 | R8 | N/A | N/A | N/A |

![Progress of in vitro selection - OG](images/Progress_of_in_vitro_selection.png)

## Progress of in vitro selection - Update (30 August 2024)

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

### Forward Paird Reads - Filtered by Structure

| Stat | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | D1 | D2 | D3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Total | 1485536 | 1533916 | 1649680 | 1436328 | 1937410 | 2336945 | 1229247 | 1756169 | 1226539 | 1090909 | 1656088 |
| Quality | 1125487 | 1177286 | 1267375 | 1107470 | 1482964 | 1800051 | 945991 | 1348014 | 913723 | 812516 | 1241307 |
| Unique | 836538 | 869219 | 836906 | 726372 | 799992 | 378395 | 38617 | 28478 | 705348 | 621216 | 830676 |
| Diversity% | 74 | 74 | 66 | 66 | 54 | 21 | 4 | 2 | 77 | 77 | 67 |
| Families | N/A | N/A | N/A | N/A | N/A | R6 | R7 | R8 | N/A | N/A | N/A |
