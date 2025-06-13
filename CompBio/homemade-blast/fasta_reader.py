# -*- coding: utf-8 -*-
"""
FASTA file reader for protein sequences
"""

def read_sequences(filename):
    """
    Read sequences from a FASTA file
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        dict: Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    # Extract sequence ID (first word after >)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    # Accumulate sequence lines
                    current_seq.append(line.upper())
            
            # Save last sequence
            if current_id:
                sequences[current_id] = ''.join(current_seq)
                
    except Exception as e:
        print("Error reading FASTA file {}: {}".format(filename, e))
        return {}
        
    return sequences

def validate_protein_sequence(sequence):
    """
    Validate if sequence contains only valid amino acid codes
    
    Args:
        sequence (str): Protein sequence
        
    Returns:
        bool: True if valid protein sequence
    """
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    return all(aa in valid_aa for aa in sequence)

def get_sequence_info(sequences):
    """
    Get basic information about loaded sequences
    
    Args:
        sequences (dict): Dictionary of sequence ID -> sequence
        
    Returns:
        dict: Basic statistics
    """
    if not sequences:
        return {}
    
    lengths = [len(seq) for seq in sequences.values()]
    
    return {
        'total_sequences': len(sequences),
        'min_length': min(lengths),
        'max_length': max(lengths),
        'avg_length': sum(lengths) / len(lengths),
        'total_residues': sum(lengths)
    } 