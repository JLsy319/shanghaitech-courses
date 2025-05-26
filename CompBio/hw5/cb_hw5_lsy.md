# Computational Biology Homework 5

## Q1

### 1.1

Needleman-Wunsch Algorithm performs global alignment. The goal is to find the optimal alignment that spans the entire length of both sequences. It is suitable for two sequences that are expected to be generally similar and of roughly same length.  
In the dynamic programming matrix, scores can be negative. The algorithm typically initializes the first row and column with gap penalties that accumulate. The optimal score is found at the bottom-right cell of the matrix.

Smith-Waterman Algorithm performs local alignment. The goal is to identify the pair of segments that have the highest similarity, regardless of the rest. It is suitable for conserved domains or motifs within divergent sequences.  
The cell score is set to zero if it would be negative. Thus, it starts from the cell with the highest score anywhere in the matrix and traces back until a cell with a score of zero is encountered.

### 1.2

Time complexity: $O(mn)$, where $m, n$ is the length of the first and the second sequence respectively.

### 1.3

Space complexity: $O(mn)$

## Q2

1. Compile a database of aligned blocks: using a database of ungapped sequence alignments representing conserved blocks from related proteins.

2. Cluster sequences based on $r$ value: $r$ value stands for the percentage that sequences are more identical to within each block. Specifically, BLOSUM62 means that all sequences within a block, sharing 62\% or more identity, are clustered.

3. The frequency of each possible amino acid pair is counted by comparing the sequences between the clusters within each block.

4. Calculate observed probabilities $q_{ij} = \frac{f_{ij}}{\sum_{k \leq l}f_{kl}}$

5. Calculate expected probabilities, given that $p_i = q_{ii} + \sum_{j \neq i} \frac{q_{ij}}{2}$

$$
e_{ij}=
\begin{cases}
i = j& p_i^2\\
i \neq j& 2p_i p_j
\end{cases}
$$

6. Calculate Log-odds ratios $S_{ij} = log_2 (\frac{q_{ij}}{e_{ij}})$ and scale it

## Q3

Tryptophan

The scores along the diagonal of a BLOSUM matrix represent the likelihood of an amino acid aligning with itself during evolution, where tryptophan has the highest score.

## Q4

BLOSUM45 matrix is better.

Distant homologs are proteins that share a common ancestor but have diverged significantly over a long evolutionary period. They typically exhibit a lower percentage of sequence identity. A BLOSUM matrix with a lower number is tuned to detect the patterns of amino acid conservation and substitution characteristic of longer evolutionary divergence.
