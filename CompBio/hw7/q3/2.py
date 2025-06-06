import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
import matplotlib.pyplot as plt

# Parameters
PDB_IDS = ["5I34", "6ZXQ", "1IWE", "1P9B", "6FM0", "7ODX", "7VF6"]
PDB_DIRECTORY = "./hw7/q3/pdb_files/"
OUTPUT_DIRECTORY = "./hw7/q3/contact_maps_output/"

CA_THRESHOLD = 8.0  # C-alpha Threshold
CB_THRESHOLD = 7.0  # C-beta Threshold

os.makedirs(PDB_DIRECTORY, exist_ok=True)
os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)

print("-" * 30)

class ChainASelect(Select):
    """只选择A链"""
    def accept_chain(self, chain):
        return chain.id == 'A'

    def accept_residue(self, residue):
        # 排除 HETATM 残基 (如水分子、配体)，除非它们是标准氨基酸的一部分
        # " " 表示是标准残基, "H_" 开头表示 HETATM, "W" 通常是水
        return residue.id[0] == ' ' # 只选择标准残基

    def accept_atom(self, atom):
        return True


def get_coords(chain, atom_type="CA"):
    """从链中获取指定类型原子的坐标，Glycine的CB特殊处理"""
    coords = []
    residue_ids = []
    for residue in chain:
        if residue.id[0] != ' ' or residue.get_resname() > "ZZZ":
            continue
        atom_coord = np.array([np.nan, np.nan, np.nan])
        if atom_type == "CB" and residue.get_resname() == "GLY":
            if "CA" in residue:
                atom_coord = residue["CA"].get_coord()
        elif atom_type in residue:
            atom_coord = residue[atom_type].get_coord()
        coords.append(atom_coord)
        residue_ids.append(f"{residue.get_resname()}{residue.id[1]}")

    return np.array(coords), residue_ids


def calculate_contact_map(coords, threshold, residue_ids, pdb_id, atom_type):
    """计算并可视化接触图"""
    num_residues = coords.shape[0]
    if num_residues == 0:
        print(f"Warning: No {atom_type} coordinates found for {pdb_id} Chain A.")
        return None

    valid_indices = ~np.isnan(coords).any(axis=1)
    coords_filtered = coords[valid_indices]
    residue_ids_filtered = [residue_ids[i] for i in range(num_residues) if valid_indices[i]]
    
    num_residues_filtered = coords_filtered.shape[0]
    if num_residues_filtered == 0:
        print(f"Warning: No valid {atom_type} coordinates found after filtering.")
        return None

    dist_matrix = np.sqrt(np.sum((coords_filtered[:, np.newaxis, :] - coords_filtered[np.newaxis, :, :])**2, axis=2))
    contact_map = (dist_matrix < threshold).astype(int)
    
    # 可选: 去除对角线和近邻接触
    # for i in range(contact_map.shape[0]):
    #     for k in range(-2, 3): # Remove contact between i and i±k (k=0,1,2)
    #         if 0 <= i + k < contact_map.shape[0]:
    #             contact_map[i, i+k] = 0
    #         if 0 <= i - k < contact_map.shape[0]:
    #             contact_map[i, i-k] = 0
    
    plt.figure(figsize=(8, 7))
    plt.imshow(contact_map, cmap='binary', origin='lower') # origin='lower' 使矩阵左下角为(0,0)
    plt.title(f"Contact Map: {pdb_id} Chain A - {atom_type} (Threshold: {threshold}Å)")
    plt.xlabel("Residue Index")
    plt.ylabel("Residue Index")
    
    output_filename = os.path.join(OUTPUT_DIRECTORY, f"{pdb_id}_chainA_{atom_type}_contact_map.png")
    plt.savefig(output_filename)
    plt.close()
    print(f"Contact Map Saved: {output_filename}")
    
    return contact_map

parser = PDBParser(QUIET=True)

for pdb_id in PDB_IDS:
    pdb_filename = os.path.join(PDB_DIRECTORY, f"{pdb_id}.pdb")
    
    if not os.path.exists(pdb_filename):
        print(f"{pdb_filename} not found")
        continue
        
    print(f"\nProcessing {pdb_id}...")
    
    try:
        structure = parser.get_structure(pdb_id, pdb_filename)
    except Exception as e:
        print(f"Error: {e}")
        continue

    if len(structure) > 1:
        print(f"Warning: {pdb_id} has multiple models, using the first one.")
    model = structure[0] 

    chain_a_exists = 'A' in [chain.id for chain in model]
    
    if not chain_a_exists:
        print(f"Warning: Chain A not found in {pdb_id}. Skipping...")
        continue

    chain_a = model['A']
    
    # --- C-alpha Contact Map ---
    print(f"Calculating {pdb_id} Chain A C-alpha contact map...")
    ca_coords, ca_residue_ids = get_coords(chain_a, "CA")
    if ca_coords.shape[0] > 0 :
        ca_contact_map = calculate_contact_map(ca_coords, CA_THRESHOLD, ca_residue_ids, pdb_id, "CA")
        if ca_contact_map is not None:
            np.save(os.path.join(OUTPUT_DIRECTORY, f"{pdb_id}_chainA_CA_contact_map.npy"), ca_contact_map)
    else:
        print(f"Warning: C-alpha coordinates not found for {pdb_id} Chain A.")

    # --- C-beta Contact Map ---
    print(f"Calculating {pdb_id} Chain A C-beta contact map...")
    cb_coords, cb_residue_ids = get_coords(chain_a, "CB")
    if cb_coords.shape[0] > 0:
        cb_contact_map = calculate_contact_map(cb_coords, CB_THRESHOLD, cb_residue_ids, pdb_id, "CB")
        if cb_contact_map is not None:
            np.save(os.path.join(OUTPUT_DIRECTORY, f"{pdb_id}_chainA_CB_contact_map.npy"), cb_contact_map)
    else:
        print(f"Warning: C-beta coordinates not found for {pdb_id} Chain A.")

    # --- Compare C-alpha and C-beta Contact Maps ---
    print("Comparing C-alpha and C-beta contact maps...")
    if ca_contact_map is not None and cb_contact_map is not None:
        if ca_contact_map.shape == cb_contact_map.shape:
            diff_map = np.abs(ca_contact_map - cb_contact_map)
            common_contacts = np.logical_and(ca_contact_map, cb_contact_map).sum()
            ca_only_contacts = (ca_contact_map - np.logical_and(ca_contact_map, cb_contact_map)).sum()
            cb_only_contacts = (cb_contact_map - np.logical_and(ca_contact_map, cb_contact_map)).sum()
            print(f"{pdb_id}: common contacts: {common_contacts}, "
                  f"C-alpha contacts: {ca_only_contacts}, "
                  f"C-beta contacts: {cb_only_contacts}")
        else:
            print(f"Warning: C-alpha and C-beta contact maps for {pdb_id} have different shapes, cannot compare.")
            print(f"CA shape: {ca_contact_map.shape}, CB shape: {cb_contact_map.shape}")

print("\nDone processing all PDB files.")