import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from pathlib import Path

def load_exclusion_ligands(filepath: Path) -> set:
    with open(filepath, 'r') as f:
        return {line.strip() for line in f if line.strip()}

def load_and_process_pdb_data(json_path: Path, exclusion_ligands: set) -> list:
    with open(json_path, 'r') as f:
        raw_data = json.load(f)

    processed_entries = []
    for item in raw_data:
        data = item.get('data', {})
        if not data:
            continue

        pdb_id = data.get('rcsb_id')
        entry_info = data.get('rcsb_entry_info')
        resolution_list = (entry_info or {}).get('resolution_combined')
        resolution = resolution_list[0] if isinstance(resolution_list, list) and resolution_list else None
        
        # Safely extract R-free value
        refine_list = data.get('refine')
        r_free = (refine_list[0] or {}).get('ls_R_factor_R_free') if isinstance(refine_list, list) and refine_list else None
        
        authors = data.get('audit_author', [])
        pi_name = authors[-1]['name'] if authors else None
        
        # Safely get Taxon ID
        polymer_entities_list = data.get('polymer_entities')
        taxon_id = None
        if isinstance(polymer_entities_list, list) and polymer_entities_list:
            first_entity = (polymer_entities_list[0] or {})
            source_organisms = first_entity.get('rcsb_entity_source_organism')
            if isinstance(source_organisms, list) and source_organisms:
                taxon_id = (source_organisms[-1] or {}).get('ncbi_taxonomy_id')

        # Determine if the structure is a complex with a drug-like molecule
        is_drug_complex = False
        nonpolymer_entities = data.get('nonpolymer_entities')
        if nonpolymer_entities:
            for ligand in nonpolymer_entities:
                ligand_id = (ligand.get('nonpolymer_comp', {}).get('chem_comp', {}) or {}).get('id')
                if ligand_id and ligand_id not in exclusion_ligands:
                    is_drug_complex = True
                    break
        
        processed_entries.append({
            "pdb_id": pdb_id,
            "resolution": resolution,
            "r_free": r_free,
            "pi_name": pi_name,
            "taxon_id": taxon_id,
            "is_drug_complex": is_drug_complex,
        })
        
    return processed_entries

def calculate_and_print_statistics(data: list):
    total_structures = len(data)
    drug_complex_count = sum(1 for entry in data if entry['is_drug_complex'])
    unique_taxa = {entry['taxon_id'] for entry in data if entry['taxon_id']}
    pi_counts = Counter(entry['pi_name'] for entry in data if entry['pi_name'])

    print(f"Analyzed a total of {total_structures} PDB entries.")
    print("\n--- Key Findings ---")
    print(f"  - Structures with Drug-like Ligands: {drug_complex_count} ({drug_complex_count/total_structures:.2%})")
    print(f"  - Unique Species Represented:        {len(unique_taxa)}")
    print(f"  - Unique Research Groups (PIs):      {len(pi_counts)}")

    print("\n--- Top 5 Most Prolific Research Groups ---")
    print(f"  {'Rank':<5} {'Principal Investigator':<30} {'Count':<10} {'Percentage'}")
    print(f"  {'-'*4:<5} {'-'*22:<30} {'-'*5:<10} {'-'*10}")
    for i, (pi, count) in enumerate(pi_counts.most_common(5), 1):
        percentage = (count / total_structures) * 100
        print(f"  {i:<5}. {pi:<30} | {count:<9} ({percentage:.1f}%)")
    print("="*60 + "\n")

def create_visualizations(data: list, output_dir: Path):
    plt.style.use('fivethirtyeight')

    resolutions = [entry['resolution'] for entry in data if entry['resolution'] is not None]
    if resolutions:
        mean_res = np.mean(resolutions)
        median_res = np.median(resolutions)

        plt.figure(figsize=(10, 6))
        plt.hist(resolutions, bins=np.arange(min(resolutions), max(resolutions) + 0.2, 0.2), 
                color='#6495ED', alpha=0.75, edgecolor='grey')
        plt.axvline(mean_res, color='red', linestyle='dashed', linewidth=2, label=f'Mean: {mean_res:.2f} Å')
        plt.axvline(median_res, color='green', linestyle='dashed', linewidth=2, label=f'Median: {median_res:.2f} Å')
        plt.legend()
        plt.title('Resolution Distribution of Mpro Structures', fontsize=18, pad=20)
        plt.xlabel('Resolution (Å)', fontsize=14)
        plt.ylabel('Structure Count', fontsize=14)
        plt.savefig(output_dir / 'Mpro_Resolution_Analysis.png', dpi=300, bbox_inches='tight')
        plt.close()

    r_free_values = [entry['r_free'] for entry in data if entry['r_free'] is not None]
    if r_free_values:
        mean_rfree = np.mean(r_free_values)
        median_rfree = np.median(r_free_values)

        plt.figure(figsize=(10, 6))
        plt.hist(r_free_values, bins=np.arange(min(r_free_values), max(r_free_values) + 0.02, 0.02), 
                color='#FF6347', alpha=0.75, edgecolor='grey')
        plt.axvline(mean_rfree, color='blue', linestyle='dashed', linewidth=2, label=f'Mean: {mean_rfree:.3f}')
        plt.axvline(median_rfree, color='purple', linestyle='dashed', linewidth=2, label=f'Median: {median_rfree:.3f}')
        plt.legend()
        plt.title('R-free Value Distribution of Mpro Structures', fontsize=18, pad=20)
        plt.xlabel('R-free Value', fontsize=14)
        plt.ylabel('Structure Count', fontsize=14)
        plt.savefig(output_dir / 'Mpro_Rfree_Analysis.png', dpi=300, bbox_inches='tight')
        plt.close()

def main():
    script_dir = Path(__file__).parent
    ligand_list_path = script_dir / 'artlist.txt'
    # Download the json file via the link in the report, as the file is too large to be uploaded to GitHub
    json_data_path = script_dir / 'rcsb_pdb_report.json'
    
    exclusion_ligands = load_exclusion_ligands(ligand_list_path)
    processed_data = load_and_process_pdb_data(json_data_path, exclusion_ligands)
    calculate_and_print_statistics(processed_data)
    create_visualizations(processed_data, output_dir=script_dir)

if __name__ == "__main__":
    main()