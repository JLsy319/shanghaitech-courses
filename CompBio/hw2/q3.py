import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit import DataStructs

def binary_string_to_bitvect(binary_string):
    """Converts a binary string (e.g., '01101') to an RDKit ExplicitBitVect."""
    if not isinstance(binary_string, str) or not all(c in '01' for c in binary_string):
        return None # Return None if input is not a valid binary string
    bv = DataStructs.ExplicitBitVect(len(binary_string))
    for i, bit in enumerate(binary_string):
        if bit == '1':
            bv.SetBit(i)
    return bv

def cluster_ligands(df, fingerprint_column_name, id_column_name, cutoff):
    """
    Groups ligands based on Tanimoto similarity using a leader-based algorithm.

    Args:
        df (pd.DataFrame): DataFrame containing ligand data and fingerprints.
        fingerprint_column_name (str): Name of the column with binary fingerprint strings.
        id_column_name (str): Name of the column with unique ligand identifiers.
        cutoff (float): Tanimoto similarity cutoff for grouping.

    Returns:
        list: A list of group assignments for each ligand in the DataFrame.
    """
    num_ligands = len(df)
    if num_ligands == 0:
        return []

    fingerprints = []
    valid_indices = []

    for i, fp_str in enumerate(df[fingerprint_column_name]):
        bv = binary_string_to_bitvect(fp_str)
        if bv:
            fingerprints.append(bv)
            valid_indices.append(i)
        else:
            print(f"Warning: Invalid or missing fingerprint for ligand at index {i} ({df.loc[i, id_column_name]}). Skipping for clustering.")

    if not fingerprints: 
        return ["No_Valid_FP"] * num_ligands

    # Map valid_indices back to original dataframe indices for group assignments
    # Initialize group_assignments for all ligands, including those with invalid FPs
    group_assignments = ["Unclustered_Invalid_FP"] * num_ligands
    
    # Temporary assignments for valid fingerprints
    temp_group_assignments = [-1] * len(fingerprints) # -1 means unassigned for valid FPs
    
    cluster_id_counter = 0
    
    for i in range(len(fingerprints)):
        if temp_group_assignments[i] == -1:
            temp_group_assignments[i] = cluster_id_counter
            leader_fp = fingerprints[i]
            
            # Compare this leader with all other unassigned ligands
            for j in range(i + 1, len(fingerprints)):
                if temp_group_assignments[j] == -1: # If ligand j is not yet assigned
                    similarity = DataStructs.TanimotoSimilarity(leader_fp, fingerprints[j])
                    if similarity >= cutoff:
                        temp_group_assignments[j] = cluster_id_counter
            cluster_id_counter += 1

    # Map temporary assignments back to the original dataframe structure
    for i, temp_idx in enumerate(valid_indices):
        group_assignments[temp_idx] = f"Group_{temp_group_assignments[i]}"
        
    return group_assignments

if __name__ == "__main__":
    input_excel_path = "./hw2/q2_output.xlsx"
    output_excel_path = "./hw2/q3.xlsx"
    
    # Define the fingerprint type to use for clustering
    # We used 'ECFP4_Binary' in the previous step.
    # Ensure this column name matches your Excel file.
    fp_column_to_use = 'ECFP4_Binary'
    ligand_id_column = 'Ligand-ID' # Or 'SMILES' or index, if preferred

    try:
        df_ligands = pd.read_excel(input_excel_path)
        print(f"Successfully read {len(df_ligands)} ligands from {input_excel_path}\n")
    except FileNotFoundError:
        print(f"Error: Input Excel file not found at {input_excel_path}")
        exit()
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        exit()

    if fp_column_to_use not in df_ligands.columns:
        print(f"Error: Fingerprint column '{fp_column_to_use}' not found in the Excel file.")
        print(f"Available columns are: {df_ligands.columns.tolist()}")
        exit()
    if ligand_id_column not in df_ligands.columns:
        print(f"Error: ID column '{ligand_id_column}' not found.")
        exit()

    # --- Clustering with Tanimoto cutoff = 0.5 ---
    print("Clustering with Tanimoto cutoff = 0.5...")
    df_ligands['Group_T0.5'] = cluster_ligands(df_ligands, fp_column_to_use, ligand_id_column, 0.5)
    print("Clustering with Tanimoto cutoff = 0.5 finished.\n")

    # --- Clustering with Tanimoto cutoff = 0.3 ---
    print("Clustering with Tanimoto cutoff = 0.3...")
    df_ligands['Group_T0.3'] = cluster_ligands(df_ligands, fp_column_to_use, ligand_id_column, 0.3)
    print("Clustering with Tanimoto cutoff = 0.3 finished.\n")

    # Display results
    print("--- Ligand Clustering Results ---")
    if 'SMILES' in df_ligands.columns:
        print(df_ligands[[ligand_id_column, 'SMILES', 'Group_T0.5', 'Group_T0.3']].head())
    else:
        print(df_ligands[[ligand_id_column, 'Group_T0.5', 'Group_T0.3']].head())

    try:
        df_ligands.to_excel(output_excel_path, index=False)
        print(f"\nSuccessfully saved clustered ligand data to {output_excel_path}")
    except Exception as e:
        print(f"\nError saving updated data to Excel: {e}")

    print("\n--- Group Counts (Tanimoto Cutoff = 0.5) ---")
    print(df_ligands['Group_T0.5'].value_counts().sort_index())

    print("\n--- Group Counts (Tanimoto Cutoff = 0.3) ---")
    print(df_ligands['Group_T0.3'].value_counts().sort_index())