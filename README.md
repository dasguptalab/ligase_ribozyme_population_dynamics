# RNA Selection Amplification

Project for analyzing RNA sequences from in vitro selection-amplification to isolate RNA ligase ribozyme.

### Data Analysis Workflow Steps

| Step | Description | Script | Directory |
| -- | -------- | ---- | --- |
| 1 | Merge the paired-end reads for each round (maximum overlap of 81-bp) | 01\_merge.sh | analysis |
| 2 | Filter reads by quality (Q>30) and remove adapter sequences | 02\_trim.sh | analysis |
| 3 | Filter reads to keep only those that contain the constant 8-bp region flanking the 40 nt variable region | 03\_filter.sh | analysis |
| 4 | Trim reads to retain only the 40-nt variable region | 04\_clean.sh | analysis |
| 5 | Combine the merged read files with the unmerged forward read files to retain the most data possible before quality control | 05\_combine.sh | analysis |
| 6 | De-replicate sequences to identify unique sequences and their read counts | 06\_format.sh | analysis |
| 7 | Filter to keep unique sequences that have >2 reads | 06\_format.sh | analysis |
| 8 | Order the unique sequences according to their read counts | 06\_format.sh | analysis |
| 9 | Cluster round 8 unique sequences according to sequence similarity (sequence identity threshold of 0.9) | 07\_cluster.sh | analysis |
| 10 | Determine the ten most abundant sequence families in round 8 and identify their peak sequences. Ensure at least 10% dissimilarity between peaks | 08\_summarize.sh | analysis |
| 11 | Identify unique sequences in each round that have at least 90% identity to any of the round 8 peak sequences | 11\_identification.sh | analysis |
| 12 | Quantify the abundance of the unique sequences in each round | 09\_quantification.sh | analysis |
| 13 | Determine the ten most abundant unique sequences in each round and determine their similarity to the round 8 families and peak sequences | 12\_top10\_identification.sh | analysis |
| 14 | Calculate the frequency of occurrence for each nucleotide in the most abundant sequence family (Family 1) | 04\_F4A\_family\_base\_conservation.R | tables\_and\_figures |
| 15 | Identify regions in the unique sequences (with >2 reads) from each round that are complementary to the substrate | 13b\_conservation.sh | analysis |
| 16 | Identify regions in the peak sequences that are complementary to the substrate in the ten most abundant families | 13d\_conservation\_families.sh | analysis |
| 17 | Identify regions in the ten most abundant sequences in each round that are complementary to the substrate | 13c\_conservation\_top10.sh | analysis |
