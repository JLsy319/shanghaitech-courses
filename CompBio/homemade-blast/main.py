#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Custom BLAST - Main Program
A fast protein sequence alignment tool with parallel processing and standard reporting

Usage:
    python main.py <database.fasta> <query.fasta> [options]

Options:
    --gapped      Perform gapped alignment (default: ungapped)
    --max-hits    Maximum hits to report per query (default: 10)
    --threshold   Minimum score threshold (default: 20)
    --kmer-size   K-mer size for indexing (default: 3)
    --threads     Number of threads for parallel search (default: all available cores)
    --output      File to save BLAST report (default: print to console)

Example:
    python main.py uniprot_sprot.fasta test_seq.fasta --gapped --threads 8 --output results.txt
"""

import sys
import argparse
import os
import multiprocessing as mp
from blast_search import BlastSearchEngine
from blast_formatter import BlastFormatter
from fasta_reader import get_sequence_info, read_sequences
import time

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Custom BLAST - Fast Protein Sequence Alignment Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Workflow:
  1. Build the database index (only needs to be done once):
     python make_blast_db.py <database.fasta>

  2. Run the search against the pre-built database:
     python main.py <database.fasta> <query.fasta> [options]
        """
    )
    
    # Required arguments
    parser.add_argument('database', help='Path to the original database FASTA file (e.g., uniprot_sprot.fasta)')
    parser.add_argument('query', help='Query FASTA file')
    
    # Optional arguments
    parser.add_argument('--gapped', action='store_true', 
                       help='Perform gapped alignment (default: ungapped only)')
    parser.add_argument('--max-hits', type=int, default=10,
                       help='Maximum hits to report per query (default: 10)')
    parser.add_argument('--threshold', type=int, default=20,
                       help='Minimum score threshold for reporting (default: 20)')
    parser.add_argument('--threads', type=int, default=mp.cpu_count(),
                       help='Number of threads for parallel processing (default: all cores)')
    parser.add_argument('--output', type=str, default=None,
                       help='Output file to save BLAST report (default: print to console)')
    
    return parser.parse_args()

def print_header():
    """Print program header"""
    print("="*80)
    print("Custom BLAST v1.1 - Parallel, Optimized Protein Alignment")
    print("="*80)

def print_parameters(args, blast_engine):
    """Print search parameters"""
    print("\nSearch Parameters:")
    print("  Database file: {}.blastdb".format(args.database))
    print("  Query file: {}".format(args.query))
    print("  Alignment type: {}".format("Gapped" if args.gapped else "Ungapped"))
    print("  Max hits: {}".format(args.max_hits))
    print("  Score threshold: {}".format(args.threshold))
    print("  K-mer size: {}".format(blast_engine.database.kmer_size))
    print("  Parallel threads: {}".format(args.threads))
    if args.output:
        print("  Output file: {}".format(args.output))

def main():
    """Main program function"""
    total_start_time = time.time()

    args = parse_arguments()
    
    print_header()
    
    try:
        # Initialize search engine with parallel processing
        print("\nInitializing BLAST search engine...")
        blast_engine = BlastSearchEngine(
            args.database, 
            num_threads=args.threads
        )
        
        print_parameters(args, blast_engine)
        
        # Get database and query info for reporting
        db_stats = blast_engine.database.get_database_stats()
        query_sequences = read_sequences(args.query)
        query_info = get_sequence_info(query_sequences)
        
        # Store all parameters for formatter.
        # Note: We can't efficiently get all subject lengths anymore without loading
        # them all from the database, so we pass an empty dict. This is a reasonable
        # trade-off for the massive performance gain in loading time.
        search_params = {
            'kmer_size': blast_engine.database.kmer_size,
            'gapped': args.gapped,
            'gap_open': blast_engine.aligner.gap_open,
            'gap_extend': blast_engine.aligner.gap_extend,
            'database_size': blast_engine.database.stats.get('total_residues'),
            'num_sequences': blast_engine.database.stats.get('total_sequences'),
            'query_length': list(query_info.values())[0] if query_info else 0,
            'subject_lengths': {} # Pass empty dict, length is shown in alignment section
        }
        
        # Perform search
        results, total_time = blast_engine.search(
            query_file=args.query,
            gapped=args.gapped,
            max_targets=args.max_hits,
            score_threshold=args.threshold
        )
        
        # Format and output results
        formatter = BlastFormatter(program_name="CustomBLAST", version="1.1")
        
        if args.output:
            formatter.save_results(results, total_time, args.database, args.query, search_params, args.output)
        else:
            # Print to console
            report = formatter.format_results(results, total_time, args.database, args.query, search_params)
            print(report)

        total_end_time = time.time()
        print("\n" + "="*80)
        print("Search completed successfully!")
        print("Total wall-clock time: {:.2f} seconds".format(total_end_time - total_start_time))
        sys.exit(0)
        
    except FileNotFoundError as e:
        print("\nError: File not found - {}".format(e))
        sys.exit(1)
    except Exception as e:
        print("\nAn error occurred: {}".format(e))
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main() 