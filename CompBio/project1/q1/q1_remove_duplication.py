import pandas as pd

def process_coronavirus_data(input_tsv: str, output_excel: str):
    """
    Process coronavirus dataset:
    1. Read TSV file
    2. Remove duplicates based on Taxon ID
    3. Save Name and Taxon ID to Excel file
    
    Args:
        input_tsv: Path to input TSV file
        output_excel: Path to output Excel file
    """
    # Read TSV file
    print(f"Reading data from {input_tsv}...")
    df = pd.read_csv(input_tsv, sep='\t')
    
    # Print initial statistics
    print(f"\nInitial dataset:")
    print(f"Total records: {len(df)}")
    print(f"Unique Taxon IDs: {df['Organism Taxonomic ID'].nunique()}")
    
    # Keep only the columns we need
    df_filtered = df[['Organism Name', 'Organism Taxonomic ID']]
    
    # Remove duplicates based on Taxon ID
    # keep='first' means we keep the first occurrence of each Taxon ID
    df_unique = df_filtered.drop_duplicates(subset=['Organism Taxonomic ID'], keep='first')
    
    # Sort by Taxon ID for better readability
    df_sorted = df_unique.sort_values(by='Organism Taxonomic ID')
    
    # Reset index after sorting
    df_final = df_sorted.reset_index(drop=True)
    
    # Rename columns to match requirements
    df_final.columns = ['Name', 'Taxon ID']
    
    # Print final statistics
    print(f"\nAfter removing duplicates:")
    print(f"Total records: {len(df_final)}")
    print(f"Unique Taxon IDs: {df_final['Taxon ID'].nunique()}")
    
    # Save to Excel
    print(f"\nSaving results to {output_excel}...")
    df_final.to_excel(output_excel, index=False)
    print("Done!")

if __name__ == "__main__":
    input_file = "q1/coronavirus_ncbi_dataset.tsv"
    output_file = "q1/q1.xlsx"
    process_coronavirus_data(input_file, output_file) 