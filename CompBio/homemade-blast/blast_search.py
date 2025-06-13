# -*- coding: utf-8 -*-
"""
Main BLAST search engine with parallel processing and optimization
"""

import time
import multiprocessing as mp
from dataclasses import dataclass, field
from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from blast_database import BlastDatabase
from blast_aligner import BlastAligner
from blast_statistics import BlastStatistics
from fasta_reader import read_sequences

# Data structures for results

@dataclass
class HSP:
    """High-scoring Segment Pair data structure (mutable)"""
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    score: float
    identity: float
    subject_id: str
    # The following fields will be filled later in the pipeline
    bit_score: float = field(default=None, repr=False)
    evalue: float = field(default=None, repr=False)
    query_coverage: float = field(default=None, repr=False)

BlastResult = namedtuple('BlastResult', ['query_id', 'hsps', 'search_time'])

# --- Worker functions for ProcessPoolExecutor ---
# These functions are defined at the top level so they can be pickled and sent to child processes.

def _process_seed_worker(query_seq, subject_seq, query_pos, subject_pos,
                         subject_id, gapped, score_threshold, aligner):
    """
    Worker function to process a single seed. This is called by _process_target_seeds_worker.
    It performs the alignment and returns an HSP object if the score is above the threshold.
    """
    # Perform ungapped extension
    q_start, q_end, s_start, s_end, score = aligner.ungapped_extension(
        query_seq, subject_seq, query_pos, subject_pos
    )

    # If gapped alignment requested and ungapped score is promising
    if gapped and score > 15:
        extend = 10
        gapped_q_start = max(0, q_start - extend)
        gapped_q_end = min(len(query_seq), q_end + extend)
        gapped_s_start = max(0, s_start - extend)
        gapped_s_end = min(len(subject_seq), s_end + extend)

        # Perform fast gapped alignment
        q_start, q_end, s_start, s_end, score = aligner.gapped_alignment_fast(
            query_seq, subject_seq,
            gapped_q_start, gapped_q_end,
            gapped_s_start, gapped_s_end
        )

    if score < score_threshold:
        return None

    identity = aligner.calculate_identity(
        query_seq, subject_seq, q_start, q_end, s_start, s_end
    )
    return HSP(q_start, q_end, s_start, s_end, score, identity, subject_id)

def _process_target_seeds_worker(query_seq, subject_seq, seed_positions, subject_id,
                                 gapped, score_threshold, aligner):
    """
    Top-level worker function for ProcessPoolExecutor.
    Processes all seeds for a single target sequence and returns a list of HSPs.
    """
    hsps = []
    for query_pos, subject_pos in seed_positions:
        hsp = _process_seed_worker(query_seq, subject_seq, query_pos, subject_pos,
                                   subject_id, gapped, score_threshold, aligner)
        if hsp:
            hsps.append(hsp)
    return hsps

class BlastSearchEngine:
    """Main BLAST search engine with parallel processing"""
    
    def __init__(self, database_file, num_threads=None):
        """
        Initialize search engine with a pre-built database
        
        Args:
            database_file (str): Path to the original FASTA database file
            num_threads (int): Number of threads for parallel processing
        """
        self.database = BlastDatabase(database_file)
        self.aligner = BlastAligner()
        self.database_file = database_file
        
        # Set number of threads for the ProcessPoolExecutor
        if num_threads is None:
            self.num_threads = mp.cpu_count()
        else:
            self.num_threads = min(num_threads, mp.cpu_count())
        
        print("Using {} threads for parallel processing".format(self.num_threads))
        
        # Initialize statistics calculator
        db_stats = self.database.get_database_stats()
        self.statistics = BlastStatistics(
            database_size=db_stats.get('total_residues', 1000000)
        )
        
    def search(self, query_file, gapped=True, max_targets=10, score_threshold=20):
        """
        Perform BLAST search
        
        Args:
            query_file (str): Path to query FASTA file
            gapped (bool): Whether to perform gapped alignment
            max_targets (int): Maximum number of hits to return per query
            score_threshold (int): Minimum score threshold for reporting hits
            
        Returns:
            tuple: (results_list, total_time)
        """
        print("\n=== Starting BLAST Search ===")
        print("Query file: {}".format(query_file))
        print("Gapped alignment: {}".format(gapped))
        print("Max targets: {}".format(max_targets))
        print("Score threshold: {}".format(score_threshold))
        
        total_start_time = time.time()
        
        # Load query sequences
        query_sequences = read_sequences(query_file)
        if not query_sequences:
            print("Error: No valid sequences found in query file")
            return [], 0
        
        print("Loaded {} query sequence(s)".format(len(query_sequences)))
        
        results = []
        
        for query_id, query_seq in query_sequences.items():
            print("\nProcessing query: {}".format(query_id))
            query_start_time = time.time()
            
            # Update statistics with query length
            self.statistics.query_length = len(query_seq)
            
            # Find all matching targets
            target_seeds = self.database.search_targets_by_seeds(query_seq)
            print("Found {} potential target sequences".format(len(target_seeds)))
            
            # Process target sequences in parallel
            hsps = self._process_targets_parallel(
                query_seq, target_seeds, gapped, score_threshold
            )
            
            # Add statistical information to HSPs
            enhanced_hsps = []
            for hsp in hsps:
                # Calculate additional statistics
                bit_score = self.statistics.calculate_bit_score(hsp.score)
                evalue = self.statistics.calculate_evalue(hsp.score, len(query_seq))
                query_coverage = self.statistics.calculate_query_coverage(
                    hsp.query_start, hsp.query_end, len(query_seq)
                )
                
                # Create enhanced HSP with additional info
                enhanced_hsp = HSP(
                    hsp.query_start, hsp.query_end, hsp.subject_start, hsp.subject_end,
                    hsp.score, hsp.identity, hsp.subject_id
                )
                # Store additional stats as attributes
                enhanced_hsp.bit_score = bit_score
                enhanced_hsp.evalue = evalue
                enhanced_hsp.query_coverage = query_coverage
                enhanced_hsps.append(enhanced_hsp)
            
            # Sort by E-value and take top results
            enhanced_hsps.sort(key=lambda h: h.evalue)
            enhanced_hsps = enhanced_hsps[:max_targets]
            
            query_time = time.time() - query_start_time
            results.append(BlastResult(query_id, enhanced_hsps, query_time))
            
            print("Query {} completed in {:.2f}s, found {} HSPs".format(
                query_id, query_time, len(enhanced_hsps)))
        
        total_time = time.time() - total_start_time
        print("\n=== BLAST Search Completed ===")
        print("Total time: {:.2f}s".format(total_time))
        
        return results, total_time
    
    def _process_targets_parallel(self, query_seq, target_seeds, gapped, score_threshold):
        """
        Process target sequences in parallel using a ProcessPoolExecutor to leverage multiple CPU cores.
        """
        MAX_PENDING = max(self.num_threads * 4, 64)
        all_hsps = []
        pending = []

        # Use ProcessPoolExecutor for true multi-core CPU-bound parallelism
        with ProcessPoolExecutor(max_workers=self.num_threads) as executor:
            for subject_id, seed_positions in target_seeds.items():
                subject_seq = self.database.get_sequence(subject_id)
                if not subject_seq:
                    continue

                clustered_seeds = self.aligner.cluster_seeds(seed_positions)
                if not clustered_seeds:
                    continue

                future = executor.submit(
                    _process_target_seeds_worker, # Use the top-level worker function
                    query_seq,
                    subject_seq,
                    clustered_seeds,
                    subject_id,
                    gapped,
                    score_threshold,
                    self.aligner  # Pass the aligner instance to the worker
                )
                pending.append(future)

                if len(pending) >= MAX_PENDING:
                    for done in as_completed(pending):
                        try:
                            hsps = done.result()
                            if hsps:
                                all_hsps.extend(hsps)
                        except Exception as e:
                            print(f"An error occurred in a worker process: {e}")
                    pending.clear()

            # Process the remaining futures
            for done in as_completed(pending):
                try:
                    hsps = done.result()
                    if hsps:
                        all_hsps.extend(hsps)
                except Exception as e:
                    print(f"An error occurred in a worker process: {e}")

        return self._filter_and_sort_hsps(all_hsps)
    
    def _filter_and_sort_hsps(self, hsps):
        """
        Remove overlapping HSPs and sort by score
        
        Args:
            hsps (list): List of HSP objects
            
        Returns:
            list: Filtered and sorted HSPs
        """
        if not hsps:
            return []
        
        # Sort by score (descending)
        hsps.sort(key=lambda x: x.score, reverse=True)
        
        # Simple overlap removal (keep highest scoring for each subject)
        seen_subjects = set()
        filtered_hsps = []
        
        for hsp in hsps:
            # For now, keep only best HSP per subject sequence
            if hsp.subject_id not in seen_subjects:
                filtered_hsps.append(hsp)
                seen_subjects.add(hsp.subject_id)
        
        return filtered_hsps
    
    def print_results(self, results, total_time):
        """
        Print formatted search results
        
        Args:
            results (list): List of BlastResult objects
            total_time (float): Total search time in seconds
        """
        print("\n" + "="*80)
        print("BLAST SEARCH RESULTS")
        print("="*80)
        print("Total search time: {:.2f} seconds".format(total_time))
        
        for result in results:
            print("\nQuery: {}".format(result.query_id))
            print("Search time: {:.2f}s".format(result.search_time))
            print("Number of hits: {}".format(len(result.hsps)))
            print("-" * 60)
            
            if not result.hsps:
                print("No significant hits found.")
                continue
            
            for i, hsp in enumerate(result.hsps[:5], 1):  # Show top 5 hits
                print("Hit {}: {}".format(i, hsp.subject_id))
                print("  Score: {:.1f}".format(hsp.score))
                print("  Identity: {:.1f}%".format(hsp.identity))
                print("  Query range: {}-{}".format(hsp.query_start, hsp.query_end))
                print("  Subject range: {}-{}".format(hsp.subject_start, hsp.subject_end))
                print()
        
        # Performance summary
        if results:
            avg_time = total_time / len(results)
            print("Average time per query: {:.2f}s".format(avg_time)) 