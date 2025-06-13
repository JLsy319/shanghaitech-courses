# -*- coding: utf-8 -*-
"""
BLAST alignment algorithms: ungapped and gapped extensions with NumPy acceleration
"""

import numpy as np
from numba import njit
from scoring_matrix import get_score, BLOSUM62

@njit
def _ungapped_extension_numba(query_indices, subject_indices, query_pos, subject_pos, 
                             word_size, x_dropoff, scoring_matrix):
    """
    Numba-accelerated version of ungapped extension.
    Operates on integer-indexed sequences for max performance.
    """
    # Calculate initial seed score
    score = 0
    for i in range(word_size):
        if (query_pos + i < len(query_indices) and 
            subject_pos + i < len(subject_indices)):
            score += scoring_matrix[query_indices[query_pos + i], subject_indices[subject_pos + i]]
    
    max_score = score
    best_query_end = query_pos + word_size
    best_subject_end = subject_pos + word_size
    
    # Extend to the right
    i = word_size
    while (query_pos + i < len(query_indices) and subject_pos + i < len(subject_indices)):
        aa_score = scoring_matrix[query_indices[query_pos + i], subject_indices[subject_pos + i]]
        score += aa_score
        
        if score > max_score:
            max_score = score
            best_query_end = query_pos + i + 1
            best_subject_end = subject_pos + i + 1
        elif score < max_score - x_dropoff:
            break
        i += 1
    
    # Extend to the left
    query_start = query_pos
    subject_start = subject_pos
    score = max_score
    
    i = 1
    while (query_pos - i >= 0 and subject_pos - i >= 0):
        aa_score = scoring_matrix[query_indices[query_pos - i], subject_indices[subject_pos - i]]
        score += aa_score
        
        if score > max_score:
            max_score = score
            query_start = query_pos - i
            subject_start = subject_pos - i
        elif score < max_score - x_dropoff:
            break
        i += 1
    
    return query_start, best_query_end, subject_start, best_subject_end, max_score

@njit
def _gapped_alignment_numba(query_indices, subject_indices, scoring_matrix, gap_open, gap_extend):
    """
    Numba-accelerated Smith-Waterman gapped alignment.
    """
    m, n = len(query_indices), len(subject_indices)
    if m == 0 or n == 0:
        return 0, 0, 0, 0, 0

    dp = np.zeros((m + 1, n + 1), dtype=np.float32)
    traceback = np.zeros((m + 1, n + 1), dtype=np.int8)  # 0: stop, 1: diag, 2: up, 3: left

    max_score = 0.0
    max_i, max_j = 0, 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = dp[i-1, j-1] + scoring_matrix[query_indices[i-1], subject_indices[j-1]]
            delete_score = dp[i-1, j] + gap_open # A gap in subject
            insert_score = dp[i, j-1] + gap_open # A gap in query

            scores = np.array([0.0, match_score, delete_score, insert_score])
            best_choice = np.argmax(scores)
            dp[i, j] = scores[best_choice]
            traceback[i, j] = best_choice

            if dp[i, j] > max_score:
                max_score = dp[i, j]
                max_i, max_j = i, j

    # Traceback from the cell with the highest score
    align_i, align_j = max_i, max_j
    while traceback[align_i, align_j] != 0:
        move = traceback[align_i, align_j]
        if move == 1: # Diagonal
            align_i -= 1
            align_j -= 1
        elif move == 2: # Up
            align_i -= 1
        elif move == 3: # Left
            align_j -= 1
    
    return align_i, max_i, align_j, max_j, max_score

class BlastAligner:
    """BLAST sequence alignment engine with NumPy acceleration"""
    
    def __init__(self, gap_open=-11, gap_extend=-1, word_size=3, x_dropoff=20):
        """
        Initialize aligner with parameters
        
        Args:
            gap_open (int): Gap opening penalty
            gap_extend (int): Gap extension penalty  
            word_size (int): Word size for initial seeds
            x_dropoff (int): X-dropoff value for extension termination
        """
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.word_size = word_size
        self.x_dropoff = x_dropoff
        
        # Create amino acid to index mapping for faster scoring
        self._create_scoring_matrix()
    
    def _create_scoring_matrix(self):
        """Create numpy scoring matrix for fast lookups"""
        # Standard amino acids
        amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
        self.aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
        
        # Create 20x20 scoring matrix
        self.scoring_matrix = np.zeros((20, 20), dtype=np.int8)
        
        for i, aa1 in enumerate(amino_acids):
            for j, aa2 in enumerate(amino_acids):
                if aa1 in BLOSUM62 and aa2 in BLOSUM62[aa1]:
                    self.scoring_matrix[i, j] = BLOSUM62[aa1][aa2]
                else:
                    self.scoring_matrix[i, j] = -1
    
    def ungapped_extension(self, query, subject, query_pos, subject_pos):
        """
        Perform ungapped extension from a seed. This is a wrapper for the Numba-jitted function.
        """
        # Convert sequences to integer indices for Numba
        try:
            query_indices = np.array([self.aa_to_idx[aa] for aa in query], dtype=np.int32)
            subject_indices = np.array([self.aa_to_idx[aa] for aa in subject], dtype=np.int32)
        except KeyError:
            # Fallback for sequences with non-standard amino acids, though slower.
            return self._ungapped_extension_fallback(query, subject, query_pos, subject_pos)

        q_start, q_end, s_start, s_end, score = _ungapped_extension_numba(
            query_indices, subject_indices, query_pos, subject_pos,
            self.word_size, self.x_dropoff, self.scoring_matrix
        )
        return q_start, q_end, s_start, s_end, int(score)

    def _ungapped_extension_fallback(self, query, subject, query_pos, subject_pos):
        """
        Original pure-Python implementation for sequences with unknown characters.
        """
        # Calculate initial seed score
        score = 0
        for i in range(self.word_size):
            if (query_pos + i < len(query) and 
                subject_pos + i < len(subject)):
                score += get_score(query[query_pos + i], subject[subject_pos + i])
        
        max_score = score
        best_query_end = query_pos + self.word_size
        best_subject_end = subject_pos + self.word_size
        
        # Extend to the right
        i = self.word_size
        while (query_pos + i < len(query) and subject_pos + i < len(subject)):
            aa_score = get_score(query[query_pos + i], subject[subject_pos + i])
            score += aa_score
            
            if score > max_score:
                max_score = score
                best_query_end = query_pos + i + 1
                best_subject_end = subject_pos + i + 1
            elif score < max_score - self.x_dropoff:
                break
            
            i += 1
        
        # Extend to the left
        query_start = query_pos
        subject_start = subject_pos
        score = max_score
        
        i = 1
        while (query_pos - i >= 0 and subject_pos - i >= 0):
            aa_score = get_score(query[query_pos - i], subject[subject_pos - i])
            score += aa_score
            
            if score > max_score:
                max_score = score
                query_start = query_pos - i
                subject_start = subject_pos - i
            elif score < max_score - self.x_dropoff:
                break
            
            i += 1
        
        return query_start, best_query_end, subject_start, best_subject_end, max_score
    
    def gapped_alignment(self, query, subject, query_start, query_end, subject_start, subject_end):
        """
        Perform gapped alignment using Smith-Waterman algorithm
        
        Args:
            query (str): Query sequence
            subject (str): Subject sequence
            query_start (int): Start position in query
            query_end (int): End position in query
            subject_start (int): Start position in subject
            subject_end (int): End position in subject
            
        Returns:
            tuple: (align_query_start, align_query_end, align_subject_start, align_subject_end, score)
        """
        # Extract subsequences for alignment
        query_subseq = query[query_start:query_end]
        subject_subseq = subject[subject_start:subject_end]
        
        m, n = len(query_subseq), len(subject_subseq)
        
        if m == 0 or n == 0:
            return query_start, query_end, subject_start, subject_end, 0
        
        # Initialize DP matrix
        dp = [[0] * (n + 1) for _ in range(m + 1)]
        
        max_score = 0
        max_i, max_j = 0, 0
        
        # Fill DP matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Match/mismatch score
                match_score = dp[i-1][j-1] + get_score(query_subseq[i-1], subject_subseq[j-1])
                
                # Gap scores
                delete_score = dp[i-1][j] + self.gap_open
                insert_score = dp[i][j-1] + self.gap_open
                
                # Take maximum (including 0 for local alignment)
                dp[i][j] = max(0, match_score, delete_score, insert_score)
                
                # Track maximum score position
                if dp[i][j] > max_score:
                    max_score = dp[i][j]
                    max_i, max_j = i, j
        
        # Traceback to find alignment boundaries
        i, j = max_i, max_j
        align_query_end = query_start + i
        align_subject_end = subject_start + j
        
        # Find start of alignment
        while i > 0 and j > 0 and dp[i][j] > 0:
            i -= 1
            j -= 1
        
        align_query_start = query_start + i
        align_subject_start = subject_start + j
        
        return align_query_start, align_query_end, align_subject_start, align_subject_end, max_score
    
    def calculate_identity(self, query, subject, query_start, query_end, subject_start, subject_end):
        """
        Calculate sequence identity percentage
        
        Args:
            query (str): Query sequence
            subject (str): Subject sequence
            query_start, query_end: Query alignment boundaries
            subject_start, subject_end: Subject alignment boundaries
            
        Returns:
            float: Identity percentage
        """
        alignment_length = min(query_end - query_start, subject_end - subject_start)
        
        if alignment_length <= 0:
            return 0.0
        
        identical = 0
        for i in range(alignment_length):
            if (query_start + i < len(query) and 
                subject_start + i < len(subject) and
                query[query_start + i] == subject[subject_start + i]):
                identical += 1
        
        return (identical / alignment_length) * 100.0 if alignment_length > 0 else 0.0
    
    def gapped_alignment_fast(self, query, subject, query_start, query_end, subject_start, subject_end):
        """
        Fast gapped alignment using NumPy and a Numba-jitted core function.
        """
        query_subseq = query[query_start:query_end]
        subject_subseq = subject[subject_start:subject_end]
        
        m, n = len(query_subseq), len(subject_subseq)
        if m == 0 or n == 0:
            return query_start, query_end, subject_start, subject_end, 0

        try:
            query_indices = np.array([self.aa_to_idx[aa] for aa in query_subseq], dtype=np.int32)
            subject_indices = np.array([self.aa_to_idx[aa] for aa in subject_subseq], dtype=np.int32)
        except KeyError:
            # Fallback for sequences with non-standard amino acids
            return self.gapped_alignment(query, subject, query_start, query_end, subject_start, subject_end)
        
        # Call the Numba-accelerated function
        a_start, a_end, s_start, s_end, score = _gapped_alignment_numba(
            query_indices, subject_indices, self.scoring_matrix, self.gap_open, self.gap_extend
        )

        return (query_start + a_start, query_start + a_end,
                subject_start + s_start, subject_start + s_end,
                float(score))
    
    def cluster_seeds(self, seeds, diagonal_threshold=5):
        """
        Cluster seeds that are on the same diagonal to avoid redundant extensions
        
        Args:
            seeds (list): List of (query_pos, subject_pos) tuples
            diagonal_threshold (int): Maximum distance to consider seeds on same diagonal
            
        Returns:
            list: List of representative seeds after clustering
        """
        if not seeds:
            return []
        
        # Group seeds by diagonal (query_pos - subject_pos)
        diagonal_groups = {}
        
        for query_pos, subject_pos in seeds:
            diagonal = query_pos - subject_pos
            if diagonal not in diagonal_groups:
                diagonal_groups[diagonal] = []
            diagonal_groups[diagonal].append((query_pos, subject_pos))
        
        # Select representative seed from each diagonal group
        clustered_seeds = []
        
        for diagonal, group_seeds in diagonal_groups.items():
            if len(group_seeds) == 1:
                clustered_seeds.append(group_seeds[0])
            else:
                # Sort by query position and take the middle one as representative
                group_seeds.sort()
                mid_idx = len(group_seeds) // 2
                clustered_seeds.append(group_seeds[mid_idx])
        
        return clustered_seeds 