import requests
import csv
import os

pdb_ids = [
    "1BBZ", "2E2B", "2F4J", "2FO0", "2G1T", "2G2H", "2GQG", "2HIW", "2HYY", "2HZ0",
    "2HZI", "2O88", "3CS9", "3EG0", "3EG1", "3EG2", "3EG3", "3EGU", "3K2M", "3PYY",
    "3QRI", "3QRJ", "3QRK", "3T04", "3UE4", "3UYO", "4J9B", "4J9C", "4J9D", "4J9E",
    "4J9F", "4J9G", "4J9H", "4J9I", "4JJB", "4JJC", "4JJD", "4TWP", "4WA9", "4ZOG",
    "5DC0", "5DC4", "5DC9", "5HU9", "5MO4", "5NP2", "5OAZ", "6NPE", "6NPU", "6NPV",
    "7DT2", "7N9G", "7PVQ", "7PVR", "7PVS", "7PVV", "7PW2", "7W7X", "7W7Y", "8H7F",
    "8H7H", "8I7S", "8I7Z"
]

def fetch_all_ligands(pdb_id):
    """修正字段路径后的配体提取函数"""
    try:
        entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        entry_response = requests.get(entry_url)
        entry_response.raise_for_status()
        entry_data = entry_response.json()
        
        ligands = entry_data.get("nonpolymer_entities", [])
        all_ligands = []
        
        for ligand in ligands:
            # 关键修正：使用正确的字段路径
            container_info = ligand.get("rcsb_nonpolymer_entity_container_identifiers", {})
            nonpoly_info = ligand.get("pdbx_entity_nonpoly", {})
            chem_comp_info = ligand.get("chem_comp", {})
            
            ligand_data = {
                "PDB ID": pdb_id,
                "Ligand ID": container_info.get("comp_id", "N/A"),
                "Name": nonpoly_info.get("name", "N/A"),
                "Formula": chem_comp_info.get("formula", "N/A")
            }
            all_ligands.append(ligand_data)
        
        return all_ligands
    except Exception as e:
        print(f"处理 {pdb_id} 时发生错误: {str(e)}")
        return []

def main():
    all_results = []
    for pdb_id in pdb_ids:
        print(f"处理 {pdb_id}...")
        ligands = fetch_all_ligands(pdb_id)
        all_results.extend(ligands)
        print(f"找到 {len(ligands)} 个配体")
    
    os.makedirs("Computational-Biology/hw2", exist_ok=True)
    csv_path = "Computational-Biology/hw2/all_ligands.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["PDB ID", "Ligand ID", "Name", "Formula"])
        writer.writeheader()
        writer.writerows(all_results)
    print(f"保存完成！共找到 {len(all_results)} 个配体")

if __name__ == "__main__":
    main()