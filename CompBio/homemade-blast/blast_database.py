# -*- coding: utf-8 -*-
"""
BLAST Database loader for SQLite-based, random-access databases.
"""

import time
import pickle
import os
import sqlite3
from collections import defaultdict

class BlastDatabase:
    """
    Manages access to a pre-built BLAST database stored in SQLite and a
    flat sequence data file. This enables extremely fast, low-memory startup.
    """
    
    def __init__(self, database_path):
        """
        Initializes the database by connecting to the SQLite file and opening
        the sequence data file for random access.
        """
        db_file = f"{database_path}.sqlite"
        self.seq_dat_path = f"{database_path}.seq_dat"

        if not os.path.exists(db_file) or not os.path.exists(self.seq_dat_path):
            raise FileNotFoundError(
                f"Database files not found. Please build the database first by running:\n"
                f"  python make_blast_db.py {database_path}\n"
            )

        print(f"Connecting to database: {db_file}")
        start_time = time.time()

        self.conn = sqlite3.connect(db_file)
        self.cursor = self.conn.cursor()

        # Load metadata into memory
        self.stats = {}
        for key, value in self.cursor.execute("SELECT key, value FROM metadata"):
            self.stats[key] = value
        self.kmer_size = int(self.stats.get('kmer_size', 3))
        
        self.seq_dat_handle = open(self.seq_dat_path, 'r')
        
        load_time = time.time() - start_time
        print(f"Database connection established and metadata loaded in {load_time:.2f}s.")

    def get_sequence(self, seq_id):
        """
        Retrieves a sequence by querying the SQLite DB for its metadata
        and then reading it from the flat data file.
        """
        self.cursor.execute("SELECT offset, length FROM sequences WHERE id=?", (seq_id,))
        result = self.cursor.fetchone()
        if result is None:
            return None
        
        offset, length = result
        self.seq_dat_handle.seek(offset)
        return self.seq_dat_handle.read(length)

    def find_seeds(self, query_sequence):
        """
        Finds k-mer matches by querying the SQLite database.
        """
        seeds = []
        for i in range(len(query_sequence) - self.kmer_size + 1):
            kmer = query_sequence[i:i + self.kmer_size]
            self.cursor.execute("SELECT locations FROM kmer_index WHERE kmer=?", (kmer,))
            result = self.cursor.fetchone()
            if result:
                locations = pickle.loads(result[0])
                for seq_id, pos in locations:
                    seeds.append((i, seq_id, pos))
        return seeds

    def search_targets_by_seeds(self, query_sequence, distance_threshold=40):
        """Finds seeds to extend based on the 'two-hit' algorithm."""
        all_seeds = self.find_seeds(query_sequence)
        diagonal_hits = defaultdict(list)
        for query_pos, subject_id, subject_pos in all_seeds:
            diagonal = query_pos - subject_pos
            diagonal_hits[(subject_id, diagonal)].append(query_pos)

        seeds_to_extend = defaultdict(list)
        for (subject_id, diagonal), query_positions in diagonal_hits.items():
            if len(query_positions) < 2:
                continue
            query_positions.sort()
            for i in range(len(query_positions) - 1):
                if (query_positions[i+1] - query_positions[i]) < distance_threshold:
                    s_pos = query_positions[i+1] - diagonal
                    seeds_to_extend[subject_id].append((query_positions[i+1], s_pos))
        
        return {sid: list(set(pos)) for sid, pos in seeds_to_extend.items()}
    
    def get_database_stats(self):
        """Returns a copy of the loaded database statistics."""
        return self.stats.copy()

    def __del__(self):
        """Ensures database connection and file handles are closed."""
        if hasattr(self, 'conn'):
            self.conn.close()
        if hasattr(self, 'seq_dat_handle'):
            self.seq_dat_handle.close() 