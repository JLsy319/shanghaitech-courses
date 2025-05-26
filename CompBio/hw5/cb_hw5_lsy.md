# Computational Biology Homework 5

## Q1

### 1.1

**1. Difference between Needleman-Wunsch (NW) and Smith-Waterman (SW) algorithms:**  
- **Alignment type**:  
  - NW performs **global alignment**, aligning sequences from end-to-end.  
  - SW performs **local alignment**, identifying regions of high similarity within subsequences.  
- **Scoring and initialization**:  
  - NW allows negative scores and initializes the first row/column with gap penalties.  
  - SW sets negative scores to **0** during matrix computation (no negative values retained) and initializes the first row/column to 0.  
- **Traceback**:  
  - NW traces back from the **bottom-right corner** (end of sequences) to the top-left.  
  - SW traces back from the **highest-scoring cell** in the matrix until encountering a 0.  

### 1.2

Time complexity of Smith-Waterman: `O(mn)` (The product of the length of the two sequences)

### 1.3

**3. Space complexity of Smith-Waterman algorithm**:  
The space complexity is **O(mn)** in the standard implementation. While optimizations (e.g., storing only two rows) can reduce space to *O(min(m,n))* for matrix filling, traceback requires retaining the full matrix to identify the optimal local alignment path. Thus, the practical space complexity remains *O(mn)*.  

---

**Summary**:  
1. NW (global) vs. SW (local): Scope, scoring rules, and traceback differ.  
2. Time: *O(mn)*.  
3. Space: *O(mn)* (full matrix retained for traceback).

// TODO: Hw5