#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Builds a pre-compiled BLAST database from a FASTA file into an SQLite database.
This approach enables fast, random access with minimal memory footprint,
suitable for extremely large databases like nr.

Usage:
    python make_blast_db.py <fasta_file> [kmer_size]
"""

import sys
import time
import sqlite3
import pickle
from collections import defaultdict
from fasta_reader import read_sequences, get_sequence_info

def create_blast_db_sqlite(fasta_file, kmer_size=3):
    """
    Reads a FASTA file and builds an SQLite database containing all necessary
    indexes and data for the BLAST search.
    """
    db_path = f"{fasta_file}.sqlite"
    seq_dat_path = f"{fasta_file}.seq_dat"
    
    print(f"Reading sequences from {fasta_file}...")
    sequences = read_sequences(fasta_file)
    if not sequences:
        print("Error: No sequences found in FASTA file.")
        return

    db_stats = get_sequence_info(sequences)
    print(f"Loaded {db_stats['total_sequences']} sequences with {db_stats['total_residues']} total residues.")

    # --- Setup SQLite database and tables ---
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Drop tables if they exist to ensure a fresh build
    c.execute("DROP TABLE IF EXISTS metadata")
    c.execute("DROP TABLE IF EXISTS sequences")
    c.execute("DROP TABLE IF EXISTS kmer_index")
    
    # Create tables
    c.execute("""
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        )
    """)
    c.execute("""
        CREATE TABLE sequences (
            id TEXT PRIMARY KEY,
            offset INTEGER NOT NULL,
            length INTEGER NOT NULL
        )
    """)
    c.execute("""
        CREATE TABLE kmer_index (
            kmer TEXT NOT NULL,
            locations BLOB NOT NULL
        )
    """)
    
    # --- 1. Populate metadata and sequences tables, create flat data file ---
    print("Writing sequence data and metadata to database...")
    current_offset = 0
    seq_data_to_insert = []
    with open(seq_dat_path, 'w') as f_dat:
        for seq_id, sequence in sequences.items():
            seq_len = len(sequence)
            seq_data_to_insert.append((seq_id, current_offset, seq_len))
            f_dat.write(sequence)
            current_offset += seq_len
    
    c.executemany("INSERT INTO sequences (id, offset, length) VALUES (?, ?, ?)", seq_data_to_insert)
    c.execute("INSERT INTO metadata (key, value) VALUES (?, ?)", ('kmer_size', str(kmer_size)))
    for key, value in db_stats.items():
        c.execute("INSERT INTO metadata (key, value) VALUES (?, ?)", (key, str(value)))

    # --- 2. Build and populate the k-mer index ---
    print(f"Building k-mer index (k={kmer_size})...")
    kmer_to_locations = defaultdict(list)
    for seq_id, sequence in sequences.items():
        for i in range(len(sequence) - kmer_size + 1):
            kmer = sequence[i:i + kmer_size]
            kmer_to_locations[kmer].append((seq_id, i))

    print(f"Index contains {len(kmer_to_locations)} unique k-mers. Writing to database...")
    
    # Serialize the list of locations using pickle and batch insert
    kmer_data_to_insert = [
        (kmer, pickle.dumps(locs)) for kmer, locs in kmer_to_locations.items()
    ]
    c.executemany("INSERT INTO kmer_index (kmer, locations) VALUES (?, ?)", kmer_data_to_insert)

    # --- 3. Create an index on the kmer column for fast lookups ---
    print("Creating database index on k-mers for fast lookups...")
    c.execute("CREATE INDEX idx_kmer ON kmer_index(kmer)")

    # --- Finalize ---
    conn.commit()
    conn.close()

    print("\nSQLite database created successfully:")
    print(f"  - {db_path}")
    print(f"  - {seq_dat_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python make_blast_db.py <fasta_file> [kmer_size]")
        sys.exit(1)

    fasta_path = sys.argv[1]
    k_size = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    
    start_time = time.time()
    create_blast_db_sqlite(fasta_path, k_size)
    end_time = time.time()
    print(f"\nTotal time to build database: {end_time - start_time:.2f} seconds.")
