# -*- coding: utf-8 -*-
"""
BLAST output formatter for standard BLAST report format
"""

import datetime
from blast_statistics import BlastStatistics

class BlastFormatter:
    """Format BLAST results in standard BLAST output format"""
    
    def __init__(self, program_name="Custom BLAST", version="1.0"):
        """
        Initialize formatter
        
        Args:
            program_name (str): Name of BLAST program
            version (str): Version number
        """
        self.program_name = program_name
        self.version = version
    
    def format_results(self, results, total_time, database_file, query_file, parameters):
        """
        Format complete BLAST results
        
        Args:
            results (list): List of BlastResult objects
            total_time (float): Total search time
            database_file (str): Database file path
            query_file (str): Query file path
            parameters (dict): Search parameters
            
        Returns:
            str: Formatted BLAST output
        """
        output = []
        
        # Header
        output.append(self._format_header(database_file, query_file, parameters))
        
        for result in results:
            # Query header
            output.append(self._format_query_header(result.query_id))
            
            if not result.hsps:
                output.append("\n***** No hits found *****\n")
                continue
            
            # Hit list summary
            output.append(self._format_hit_list_summary(result.hsps, result.query_id, parameters))
            
            # Alignments section
            output.append(self._format_alignments(result.hsps, result.query_id, parameters))
        
        # Footer
        output.append(self._format_footer(parameters, total_time))
        
        return '\n'.join(output)
    
    def _format_header(self, database_file, query_file, parameters):
        """Format the header section"""
        header = []
        header.append("="*80)
        header.append("{} {}".format(self.program_name, self.version))
        header.append("="*80)
        header.append("")
        header.append("Query= {}".format(query_file))
        header.append("Database: {}".format(database_file))
        header.append("")
        header.append("Search Parameters:")
        header.append("  Word size: {}".format(parameters.get('kmer_size', 3)))
        header.append("  Gap open: {}".format(parameters.get('gap_open', -11)))
        header.append("  Gap extend: {}".format(parameters.get('gap_extend', -1)))
        header.append("  Matrix: BLOSUM62")
        header.append("  Gapped: {}".format("Yes" if parameters.get('gapped', False) else "No"))
        header.append("")
        
        return '\n'.join(header)
    
    def _format_query_header(self, query_id):
        """Format query header"""
        return "\nQuery= {}\n".format(query_id)
    
    def _format_hit_list_summary(self, hsps, query_id, parameters):
        """Format hit list summary table"""
        if not hsps:
            return ""
        
        # Calculate statistics for each HSP
        stats = BlastStatistics(
            database_size=parameters.get('database_size', 1000000),
            query_length=parameters.get('query_length', 100)
        )
        
        lines = []
        lines.append("Sequences producing significant alignments:")
        lines.append("")
        lines.append("{:<50} {:>10} {:>10} {:>10} {:>10} {:>8}".format(
            "Description", "Max Score", "Total Score", "Query Cover", "E-value", "Ident"
        ))
        lines.append("-" * 100)
        
        # Sort HSPs by E-value (ascending)
        sorted_hsps = sorted(hsps, key=lambda h: stats.calculate_evalue(h.score))
        
        for hsp in sorted_hsps:
            # Calculate statistics
            bit_score = stats.calculate_bit_score(hsp.score)
            evalue = stats.calculate_evalue(hsp.score)
            query_coverage = stats.calculate_query_coverage(
                hsp.query_start, hsp.query_end, parameters.get('query_length', 100)
            )
            
            # Format E-value
            if evalue < 1e-100:
                evalue_str = "0.0"
            elif evalue < 1e-10:
                evalue_str = "{:.1e}".format(evalue)
            else:
                evalue_str = "{:.2f}".format(evalue)
            
            # Truncate description if too long
            description = hsp.subject_id[:45] + "..." if len(hsp.subject_id) > 48 else hsp.subject_id
            
            lines.append("{:<50} {:>10.1f} {:>10.1f} {:>9.0f}% {:>10} {:>7.0f}%".format(
                description, bit_score, bit_score, query_coverage, evalue_str, hsp.identity
            ))
        
        lines.append("")
        return '\n'.join(lines)
    
    def _format_alignments(self, hsps, query_id, parameters):
        """Format detailed alignments section"""
        if not hsps:
            return ""
        
        stats = BlastStatistics(
            database_size=parameters.get('database_size', 1000000),
            query_length=parameters.get('query_length', 100)
        )
        
        lines = []
        lines.append("ALIGNMENTS")
        lines.append("="*80)
        
        for i, hsp in enumerate(hsps, 1):
            lines.append("")
            lines.append("> {}".format(hsp.subject_id))
            lines.append("Length={}".format(parameters.get('subject_lengths', {}).get(hsp.subject_id, "Unknown")))
            lines.append("")
            
            # Calculate statistics
            bit_score = stats.calculate_bit_score(hsp.score)
            evalue = stats.calculate_evalue(hsp.score)
            
            lines.append(" Score = {:.1f} bits ({}),  Expect = {:.2e}".format(
                bit_score, int(hsp.score), evalue
            ))
            lines.append(" Identities = {}/{} ({:.0f}%), Gaps = 0/{} (0%)".format(
                int(hsp.identity * (hsp.query_end - hsp.query_start) / 100),
                hsp.query_end - hsp.query_start,
                hsp.identity,
                hsp.query_end - hsp.query_start
            ))
            lines.append("")
            
            # Format alignment (simplified version)
            lines.append("Query  {:>3}  {}  {}".format(
                hsp.query_start + 1,
                "..." if hsp.query_end - hsp.query_start > 60 else "SEQUENCE",
                hsp.query_end
            ))
            lines.append("            {}".format("|||" * min(20, (hsp.query_end - hsp.query_start) // 3)))
            lines.append("Sbjct  {:>3}  {}  {}".format(
                hsp.subject_start + 1,
                "..." if hsp.subject_end - hsp.subject_start > 60 else "SEQUENCE",
                hsp.subject_end
            ))
            lines.append("")
        
        return '\n'.join(lines)
    
    def _format_footer(self, parameters, total_time):
        """Format footer with parameters and statistics"""
        footer = []
        footer.append("="*80)
        footer.append("SEARCH SUMMARY")
        footer.append("="*80)
        footer.append("")
        footer.append("Matrix: BLOSUM62")
        footer.append("Gap Penalties: Existence: {}, Extension: {}".format(
            parameters.get('gap_open', -11), parameters.get('gap_extend', -1)
        ))
        footer.append("Number of sequences in database: {}".format(
            parameters.get('num_sequences', "Unknown")
        ))
        footer.append("Number of residues in database: {}".format(
            parameters.get('database_size', "Unknown")
        ))
        footer.append("")
        footer.append("Total search time: {:.2f} seconds".format(total_time))
        footer.append("Timestamp: {}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        footer.append("")
        
        return '\n'.join(footer)
    
    def save_results(self, results, total_time, database_file, query_file, parameters, output_file):
        """
        Save formatted results to file
        
        Args:
            results (list): List of BlastResult objects
            total_time (float): Total search time
            database_file (str): Database file path
            query_file (str): Query file path
            parameters (dict): Search parameters
            output_file (str): Output file path
        """
        formatted_output = self.format_results(results, total_time, database_file, query_file, parameters)
        
        with open(output_file, 'w') as f:
            f.write(formatted_output)
        
        print("Results saved to: {}".format(output_file)) 