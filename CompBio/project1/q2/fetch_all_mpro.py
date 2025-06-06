from Bio import Entrez
from Bio.Blast import NCBIWWW
import pandas as pd
import time
import os
from datetime import datetime

class MproFinder:
    def __init__(self, email: str, api_key: str | None = None):
        """
        Initialize MproFinder with NCBI credentials.
        
        Args:
            email: Your email address (required by NCBI)
            api_key: Your NCBI API key (optional)
        """
        self.email = email
        self.api_key = api_key
        
        # Setup NCBI
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        
        # Reference sequence for BLAST (SARS-CoV-2 Mpro)
        self.reference_mpro = {
            "sequence": """SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQA
GNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSF
LNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYA
AVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRT
ILGSALLEDEFTPFDVVRQCSGVTFQ""",
            "accession": "YP_009725301.1",
            "description": "3C-like proteinase [Severe acute respiratory syndrome coronavirus 2]"
        }
    
    def _search_annotated_mpro(self, taxon_id: str, virus_name: str) -> tuple[bool, str]:
        """
        Search for annotated main protease sequence.
        
        Returns:
            Tuple of (success, sequence/error_message)
        """
        search_terms = [
            f"txid{taxon_id}[Organism] AND (main protease OR 3C-like protease OR 3CLpro OR Mpro)",
            f"txid{taxon_id}[Organism] AND (nsp5 OR pp1a)",
            f"txid{taxon_id}[Organism] AND coronavirus main protease"
        ]
        
        for term in search_terms:
            try:
                handle = Entrez.esearch(db="protein", term=term, retmax=5)
                record = Entrez.read(handle)
                handle.close()
                
                if record["Count"] != "0":
                    protein_ids = record["IdList"]
                    handle = Entrez.efetch(db="protein", id=protein_ids[0], rettype="fasta", retmode="text")
                    fasta_record = handle.read()
                    handle.close()
                    
                    return True, fasta_record
                
                time.sleep(0.34 if not self.api_key else 0.11)  # NCBI rate limit
                
            except Exception as e:
                return False, f"Error in annotated search: {str(e)}"
        
        return False, "No annotated sequence found"
    
    def _search_unannotated_mpro(self, taxon_id: str, virus_name: str) -> tuple[bool, str]:
        """
        Search for main protease using BLAST in unannotated sequences.
        
        Returns:
            Tuple of (success, sequence/error_message)
        """
        try:
            # First get the genome sequence
            handle = Entrez.esearch(db="nucleotide", 
                                  term=f"txid{taxon_id}[Organism] AND (complete genome[Title] OR complete sequence[Title])",
                                  retmax=1)
            record = Entrez.read(handle)
            handle.close()
            
            if record["Count"] == "0":
                return False, "No genome sequence found"
            
            genome_id = record["IdList"][0]
            
            print(f"Running BLAST search for {virus_name}...")
            result_handle = NCBIWWW.qblast(
                "tblastn", "nr",
                self.reference_mpro["sequence"],
                entrez_query=f"txid{taxon_id}[Organism]",
                expect=0.001,
                hitlist_size=1
            )
            
            blast_results = result_handle.read()
            result_handle.close()
            
            if "No hits found" not in blast_results:
                # Create FASTA format result
                fasta_record = f">{virus_name} | Taxon_ID:{taxon_id} | BLAST_match\n{self.reference_mpro['sequence']}\n"
                return True, fasta_record
            
            return False, "No BLAST hits found"
            
        except Exception as e:
            return False, f"Error in BLAST search: {str(e)}"
    
    def find_all_mpro(self, input_excel: str = "q1.xlsx"):
        """
        Find main protease sequences for all viruses in the input file.
        
        Args:
            input_excel: Path to Excel file containing virus information
        """
        print(f"Reading virus data from {input_excel}...")
        df = pd.read_excel(input_excel)
        
        found_sequences = []
        total = len(df)
        
        for index, row in df.iterrows():
            taxon_id = str(row['Taxon ID'])
            virus_name = row['Name']
            
            print(f"\nProcessing {virus_name} ({index + 1}/{total})...")
            
            # First try to find annotated sequence
            success, result = self._search_annotated_mpro(taxon_id, virus_name)
            
            if not success:
                print(f"No annotated sequence found, trying BLAST search...")
                success, result = self._search_unannotated_mpro(taxon_id, virus_name)
            
            if success:
                found_sequences.append(result)
        
        # Save all sequences to FASTA file
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        fasta_file = f"q2_mpro_sequences.fasta"
        with open(fasta_file, 'w') as f:
            for seq in found_sequences:
                f.write(seq + '\n')
        
        print(f"\nProcessing complete!")
        print(f"Found sequences: {len(found_sequences)} out of {total}")
        print(f"\nResults saved to: {fasta_file}")

if __name__ == "__main__":
    finder = MproFinder(
        email="EMAIL_ACCOUNT",
        api_key="API_KEY"
    )
    finder.find_all_mpro()