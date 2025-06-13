# -*- coding: utf-8 -*-
"""
BLAST statistical calculations including E-value computation
"""

import math

class BlastStatistics:
    """Statistical calculations for BLAST results"""
    
    def __init__(self, database_size, query_length=0):
        """
        Initialize with database and query information
        
        Args:
            database_size (int): Total number of residues in database
            query_length (int): Length of query sequence
        """
        # Ensure database_size is an integer
        self.database_size = int(database_size)
        self.query_length = query_length
        
        # Karlin-Altschul parameters for BLOSUM62
        self.lambda_value = 0.267
        self.k_value = 0.041
        self.h_value = 0.401
    
    def calculate_bit_score(self, raw_score):
        """
        Calculate bit score from raw score
        
        Args:
            raw_score (float): Raw alignment score
            
        Returns:
            float: Bit score
        """
        return (self.lambda_value * raw_score - math.log(self.k_value)) / math.log(2)
    
    def calculate_evalue(self, raw_score, query_length=None, database_size=None):
        """
        Calculate E-value for alignment score
        
        Args:
            raw_score (float): Raw alignment score
            query_length (int): Query sequence length (optional if set in init)
            database_size (int): Database size (optional if set in init)
            
        Returns:
            float: E-value
        """
        if query_length is None:
            query_length = self.query_length
        if database_size is None:
            database_size = self.database_size
            
        if query_length is None or database_size is None:
            return float('inf')
        
        # Calculate effective search space
        effective_search_space = query_length * database_size
        
        # Calculate E-value using Karlin-Altschul statistics
        evalue = effective_search_space * math.exp(-self.lambda_value * raw_score)
        
        return max(evalue, 1e-200)  # Prevent underflow
    
    def calculate_query_coverage(self, query_start, query_end, query_length):
        """
        Calculate query coverage percentage
        
        Args:
            query_start (int): Start position in query
            query_end (int): End position in query
            query_length (int): Total query length
            
        Returns:
            float: Coverage percentage
        """
        if query_length <= 0:
            return 0.0
        
        coverage_length = query_end - query_start
        return (coverage_length / query_length) * 100.0
    
    def get_statistical_parameters(self):
        """
        Get the statistical parameters used
        
        Returns:
            dict: Dictionary of parameters
        """
        return {
            'lambda': self.lambda_value,
            'k': self.k_value,
            'h': self.h_value,
            'database_size': self.database_size,
            'query_length': self.query_length
        } 