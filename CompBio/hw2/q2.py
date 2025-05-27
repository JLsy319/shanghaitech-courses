from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
import pandas as pd

def generate_fingerprints_from_csv(input_csv_path="q2-start.csv", output_excel_path="abl1_ligand_fingerprints_from_csv.xlsx"):
    """
    Reads ligand data from a CSV file, generates MACCS and ECFP4 fingerprints,
    and saves the results to an Excel file.

    Args:
        input_csv_path (str): Path to the input CSV file.
                              CSV should have columns: PDB_ID, Ligand_ID, SMILES
        output_excel_path (str): Path to save the output Excel file.
    """
    results = []

    try:
        df_input = pd.read_csv(input_csv_path)
        print(f"Successfully read {len(df_input)} rows from {input_csv_path}")
    except FileNotFoundError:
        print(f"Error: Input CSV file not found at {input_csv_path}")
        return
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    for index, row in df_input.iterrows():
        pdb_id = row["PDB_ID"]
        ligand_id = row["Ligand_ID"]
        smiles = row["SMILES"]

        if not isinstance(smiles, str):
            print(f"Warning: SMILES for PDB_ID {pdb_id}, Ligand_ID {ligand_id} is not a string (found {type(smiles)}: {smiles}). Skipping.")
            fp_maccs_binary = "Invalid SMILES (not a string)"
            fp_ecfp4_binary = "Invalid SMILES (not a string)"
        else:
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                print(f"Warning: Could not parse SMILES for PDB_ID {pdb_id}, Ligand_ID {ligand_id}: {smiles}. Skipping fingerprint generation.")
                fp_maccs_binary = "Error parsing SMILES"
                fp_ecfp4_binary = "Error parsing SMILES"
            else:
                # MACCS Finnerprints
                fp_maccs = MACCSkeys.GenMACCSKeys(mol)
                fp_maccs_binary = fp_maccs.ToBitString()

                # ECFP4 Fingerprints
                fp_ecfp4 = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024) # nBits = 1024 or 2048
                fp_ecfp4_binary = fp_ecfp4.ToBitString()
        results.append({
            "PDB-ID": pdb_id,
            "Ligand-ID": ligand_id,
            "SMILES": smiles,
            "MACCS_Binary": fp_maccs_binary,
            "ECFP4_Binary": fp_ecfp4_binary
        })

    df_results = pd.DataFrame(results)
    print(f"\nProcessed {len(df_results)} ligands.")
    if not df_results.empty:
        print(df_results.head()) 
    try:
        df_results.to_excel(output_excel_path, index=False)
        print(f"\nSuccessfully saved data to {output_excel_path}")
    except Exception as e:
        print(f"\nError saving to Excel: {e}")
        print("Please ensure you have 'openpyxl' installed: pip install openpyxl")

if __name__ == "__main__":
    generate_fingerprints_from_csv(input_csv_path="./hw2/q2-start.csv",
                                   output_excel_path="./hw2/q2_output.xlsx")