from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBList
import pandas as pd
import requests

def get_ligand_smiles(pdb_id, ligand_id):
    """通过PDB API获取配体SMILES"""
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data["rcsb_chem_comp_descriptor"]["smiles"]
    return None

# 输入符合条件的PDB ID列表（示例）
pdb_ids = ["4j9h", "2hyy", "3k5v"]  # 替换为实际筛选的PDB IDs

# 收集数据
data = []
for pdb_id in pdb_ids:
    # 获取结构中的配体列表（需提前筛选分子量≥250且非ATP类似物）
    structure_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    entry_data = requests.get(structure_url).json()
    ligands = entry_data["rcsb_entry_info"]["polymer_entity_count_hetatm"]
    
    for lig in ligands:
        lig_id = lig["comp_id"]
        # 获取SMILES
        smiles = get_ligand_smiles(pdb_id, lig_id)
        if smiles:
            data.append({
                "PDB-ID": pdb_id,
                "Ligand-ID": lig_id,
                "SMILES": smiles
            })

# 保存为CSV
df = pd.DataFrame(data)
df.to_csv("abl1_ligands.csv", index=False)