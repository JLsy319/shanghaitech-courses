# Computational Biology Homework 6

## Q1

The Expect value (E-value) is a parameter that describes the number of hits one can “expect” to see by chance when searching a database of a particular size. It decreases exponentially as the Score of the match increases. Essentially, the E value describes the random background noise.

Given that E-value is a number to indicate the possibility, E-value could not be larger than 1.

## Q2

1. Increasing K

   - Better specificity, less sensitivity, faster
   - Useful for finding strong matches in closely related sequences

2. Decreasing K

   - Better sensitivity, less specificity, slower
   - Useful for identifying divergent sequences in evolutionary distant sequences

3. Increasing T

   - Better stringency, faster
   - It may discard biologically relevant sequences with moderate similarity

4. Decreasing T

   - Better sensitivity, slower
   - Noisy matched sequences

## Q3

The 2-hit heuristic based on that true homologs often share multiple conserved regions in close proximity. By requiring dual seeds, BLAST minimizes computational overhead and prioritizes high-confidence regions. It keeps tha balance between speed and sensitivity by reducing extensions and capturing longer or divergent alignments.

## Q4

PSI-BLAST dynamically constructs a position-specific scoring matrix through iterative searches.  
The PSSM assigns varying weights to residues based on their conservation in aligned sequences, amplifying signals from critical motifs even in otherwise divergent regions. Each iteration refines the PSSM by incorporating newly detected homologs, enabling the algorithm to learn the recognition of remote relatives.

## Q5

Transition Probability Matrix and Emission Probability Matrix. The first one models the evolution of hidden states over time. The second one links hidden states to observable data.

## Q6

1. Evaluation Problem

   - Application: Given a profile HMM, any given path through the model will emit a sequence with an associated probability.
   - Algorithm: Forward Algorithm
   - Tool: `hmmsearch`

2. Decoding Problem

   - Application: Find the most likely sequence of states that could have produced the output sequence
   - Algorithm: Viterbi Algorithm
   - Tool: `hmmscan`

3. Learning Problem

   - Application: Build profile hmm with MSA
   - Algorithm: Baum-Welch Algorithm
   - Tool: `hmmbuild`

## Q7

The key idea of PatterHunter lies in its use of optimized spaced seeds rather than the contiguous seeds. It increases the sensitivity, considering distant homologs often retaining critical residues at non-consecutive positions, as spaced seeds target these dispersed conserved sites, whereas BLAST miss them if the exact consecutive stretch diverges.
