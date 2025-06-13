# BLAST Tool

## Usage

First, build a pre-compiled BLAST database from a `FASTA` file into an SQLite database.

```bash
python make_blast_db.py <fasta_file> [kmer_size] # default = 3
```

Input the query, it will automatically point to the SQLite datebase, though the input looks like a `FASTA` file.

```bash
python main.py <database.fasta> <query.fasta> [options]
```

Options:
- `--gapped`: Perform gapped alignment (default: ungapped)
- `--max-hits`: Maximum hits to report per query (default: 10)
- `--threshold`: Minimum score threshold (default: 20)
- `--kmer-size`: K-mer size for indexing (default: 3)
- `--threads`: Number of threads for parallel search (default: all available cores)
- `--output`: File to save BLAST report (default: print to console)

## Target

- With the ultimate goal of searching the entire nr database, and fully utilizing the multi-core computing capabilities of the CPU

- Functional implementation that does not rely on excessive external libraries

## Design 1: Scalable Database Backend

Problem: Loading the entire database (sequences + index) into memory is infeasible for large databases (e.g., nr), leading to memory exhaustion and slow startup.

Solution: Use SQLite as the database backend and separate index from data.

- `.sqlite` file: Stores K-mer index and sequence metadata.
- `.seq_dat` file: Stores raw, concatenated sequence data.

```python
# Pseudo-code for get_sequence(seq_id)
def get_sequence(self, seq_id):
    # Query metadata in SQLite
    cursor.execute("SELECT offset, length FROM sequences WHERE id=?", (seq_id,))
    offset, length = cursor.fetchone()

    # Jump directly to the exact position in the data file and only read the required sequence length.
    self.seq_data_handle.seek(offset)
    return self.seq_data_handle.read(length)
```

## Design 2: Two-Hit Algorithm

Problem: When searching large databases, the number of single K-mer matches ("seeds") is enormous, with most being noise from random collisions. Extending alignments for every seed would waste massive computational resources.

Solution: Implement a stricter "Two-Hit" seed trigger mechanism than the original BLAST.

Steps:

1. Find all single-hit K-mer seeds.

2. Group them by diagonal (diagonal = query_pos - subject_pos).

3. On each diagonal, only consider it a valid "two-hit" event if two hits are found within a short distance (e.g., < 40aa).

4. Only seeds confirmed by "two-hit" events proceed to the costly extension phase.

```python
# Pseudo-code for Two-Hit Algorithm  
def search_targets_by_seeds(query_sequence):  
    all_single_hits = find_all_kmer_matches(query_sequence)  
    diagonals = group_hits_by_diagonal(all_single_hits)  

    for diag, hits in diagonals.items():  
        if len(hits) < 2: continue  # Must have at least two hits  

        sorted_hits = sort(hits)  
        for i in range(len(sorted_hits) - 1):  
            # Check distance between adjacent hits  
            if (sorted_hits[i+1].pos - sorted_hits[i].pos) < DISTANCE_THRESHOLD:  
                # Only qualified seeds proceed to extension  
                add_to_extension_list(sorted_hits[i+1]) 
```

## Design 3: Tiered Alignment Strategy

> Similar to NCBI BLAST

### Stage 1: Ungapped Extension

- Quickly extend from "two-hit" seeds in both directions to obtain a preliminary, gap-free alignment score.
- A simple linear scan, accumulating match scores until falling below a threshold.

```python
# Pseudo-code for Ungapped Extension
score = seed_score
max_score = score
# Extend right
while score > max_score - X_DROP_OFF:
    score += get_score_at_next_pos()
    update_max_score()
# Extend left
...
```

### Stage 2: Gapped Alignment

- Perform precise, sensitive gapped alignment on high-scoring regions that pass the first stage.
- Smith-Waterman dynamic programming.

```python
# Smith-Waterman Recurrence Relation
# S[i,j] = max(
#     S[i-1,j-1] + score(q_i, s_j), // Match
#     S[i-1,j]   + gap_penalty,    // Deletion
#     S[i, j-1]  + gap_penalty,    // Insertion
#     0
# )
```

## Design 4: Statistical Significance Evaluation

> Bit Score & E-value

- Bit Score: Normalizes the raw score to be independent of the scoring matrix.

$$
S' = \frac{\lambda S - \ln K}{\ln 2}
$$

- E-value: The expected number of random matches with a score ≥ the current alignment in a database of this size. Closer to 0 means more significant.

$$
E = \text{SearchSpace} \cdot K \cdot e^{-\lambda S}
$$

where $\lambda$ and $K$ are Karlin-Altschul parameters dependent on the scoring matrix and sequence background.

```python
# from blast_statistics.py
def calculate_evalue(self, raw_score):
    # Ensure database_size is a number
    db_size = int(self.database_size) 
    query_len = int(self.query_length)

    effective_search_space = (db_size - query_len) * (query_len)
    
    evalue = effective_search_space * self.k_value * math.exp(-self.lambda_value * raw_score)
    return evalue
```

## Design 5: Speed Optimization

### JIT Compilation with Numba

Problem: Python loops are slow for compute-intensive tasks like sequence alignment.

Solution: Use Numba's @njit decorator to compile performance-critical functions (e.g., `ungapped_extension`, `gapped_alignment`) into efficient machine code.

```python
from numba import njit

@njit # Just one line of code!
def _gapped_alignment_numba(...):
    # The same double loop for DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            ... # Complex calculations
    return ...
```

### Parallel Computing with ProcessPoolExecutor

Problem: Python’s Global Interpreter Lock limits multithreading for CPU-bound tasks.

Solution: Use `ProcessPoolExecutor` to spawn independent subprocesses, each with its own interpreter and memory, bypassing the GIL for true multi-core parallelism.

***Encapsulate core alignment logic into a top-level "worker function" that can be serialized and distributed.***


```python
# Using ProcessPoolExecutor to distribute tasks
with ProcessPoolExecutor(max_workers=self.num_threads) as executor:
    for subject_id, seeds in target_seeds.items():
        # Submit tasks to worker functions
        future = executor.submit(
            _process_target_seeds_worker, 
            ..., 
            self.aligner # Pass necessary objects
        )
```

## Test

Program is tested with Swiss-Prot database.

- Download it by:
```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

- Unzip
```bash
gunzip uniprot_sprot.fasta.gz
```

- A normaltesting sequence:　[P68871](https://www.uniprot.org/uniprotkb/P68871/entry#sequences)

- A challenging testing sequence: [P02452](https://www.uniprot.org/uniprotkb/P02452/entry#sequences)