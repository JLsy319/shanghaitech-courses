import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
import matplotlib.pyplot as plt

# --- 配置参数 ---
PDB_IDS = ["5I34", "6ZXQ", "1IWE", "1P9B", "6FM0", "7ODX", "7VF6"]
PDB_DIRECTORY = "./hw7/q3/pdb_files/"  # 假设你的 PDB 文件都存放在这个子目录中
OUTPUT_DIRECTORY = "./hw7/q3/contact_maps_output/" # 输出接触图的目录

CA_THRESHOLD = 8.0  # C-alpha Threshold
CB_THRESHOLD = 7.0  # C-beta Threshold

os.makedirs(PDB_DIRECTORY, exist_ok=True)
os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)

print("-" * 30)

# --- 辅助类和函数 ---
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
    residue_ids = [] # 记录残基号，用于绘图标签
    for residue in chain:
        # 确保是标准氨基酸残基 (避免读取配体等的原子)
        if residue.id[0] != ' ' or residue.get_resname() > "ZZZ": # 标准氨基酸通常用三个大写字母
            continue

        atom_coord = np.array([np.nan, np.nan, np.nan]) # 默认值
        if atom_type == "CB" and residue.get_resname() == "GLY":
            if "CA" in residue:
                atom_coord = residue["CA"].get_coord()
        elif atom_type in residue:
            atom_coord = residue[atom_type].get_coord()
        
        coords.append(atom_coord)
        residue_ids.append(f"{residue.get_resname()}{residue.id[1]}") # 例如 ALA10

    return np.array(coords), residue_ids


def calculate_contact_map(coords, threshold, residue_ids, pdb_id, atom_type):
    """计算并可视化接触图"""
    num_residues = coords.shape[0]
    if num_residues == 0:
        print(f"警告: 在 {pdb_id} 的 A 链中没有找到 {atom_type} 原子坐标。")
        return None

    # 移除包含nan的行 (如果某些残基没有所需原子)
    valid_indices = ~np.isnan(coords).any(axis=1)
    coords_filtered = coords[valid_indices]
    residue_ids_filtered = [residue_ids[i] for i in range(num_residues) if valid_indices[i]]
    
    num_residues_filtered = coords_filtered.shape[0]
    if num_residues_filtered == 0:
        print(f"警告: 在 {pdb_id} 的 A 链中过滤后没有有效的 {atom_type} 原子坐标。")
        return None

    dist_matrix = np.sqrt(np.sum((coords_filtered[:, np.newaxis, :] - coords_filtered[np.newaxis, :, :])**2, axis=2))
    
    contact_map = (dist_matrix < threshold).astype(int)
    
    # 可选: 去除对角线和近邻接触 (例如，i 与 i, i与i+1, i与i+2)
    # for i in range(contact_map.shape[0]):
    #     for k in range(-2, 3): # 示例：去除 i 和 i±k (k=0,1,2) 之间的接触
    #         if 0 <= i + k < contact_map.shape[0]:
    #             contact_map[i, i+k] = 0
    #         if 0 <= i - k < contact_map.shape[0]: # 确保对称
    #             contact_map[i, i-k] = 0
    
    # 可视化并保存
    plt.figure(figsize=(8, 7))
    plt.imshow(contact_map, cmap='binary', origin='lower') # origin='lower' 使矩阵左下角为(0,0)
    plt.title(f"Contact Map: {pdb_id} Chain A - {atom_type} (Threshold: {threshold}Å)")
    plt.xlabel("Residue Index")
    plt.ylabel("Residue Index")
    # 如果残基数量不多，可以尝试显示残基标签
    # if num_residues_filtered < 30:
    #     plt.xticks(np.arange(num_residues_filtered), residue_ids_filtered, rotation=90, fontsize=8)
    #     plt.yticks(np.arange(num_residues_filtered), residue_ids_filtered, fontsize=8)
    
    output_filename = os.path.join(OUTPUT_DIRECTORY, f"{pdb_id}_chainA_{atom_type}_contact_map.png")
    plt.savefig(output_filename)
    plt.close()
    print(f"已保存接触图: {output_filename}")
    
    return contact_map

# --- 主处理循环 ---
parser = PDBParser(QUIET=True) # QUIET=True 避免过多警告

for pdb_id in PDB_IDS:
    pdb_filename = os.path.join(PDB_DIRECTORY, f"{pdb_id}.pdb")
    
    if not os.path.exists(pdb_filename):
        print(f"PDB 文件 {pdb_filename} 未找到，跳过...")
        continue
        
    print(f"\n正在处理 {pdb_id}...")
    
    try:
        structure = parser.get_structure(pdb_id, pdb_filename)
    except Exception as e:
        print(f"解析PDB文件 {pdb_filename} 失败: {e}")
        continue

    # 检查是否有多个模型，通常我们只用第一个模型
    if len(structure) > 1:
        print(f"警告: {pdb_id} 包含多个模型，将只使用第一个模型 (ID: {structure[0].id})。")
    model = structure[0] 

    # 提取A链 (如果PDB本身就是单链A，这一步也能正确处理)
    # 有些PDB可能本身就是纯A链，有些则需要筛选
    chain_a_exists = 'A' in [chain.id for chain in model]
    
    if not chain_a_exists:
        print(f"警告: 在 {pdb_id} 中未找到 A 链，跳过该PDB的接触图计算。")
        # 如果你知道某些PDB的链名不是'A'但你想用，需要修改这里的逻辑
        continue

    chain_a = model['A']
    
    # --- C-alpha 接触图 ---
    print(f"计算 {pdb_id} Chain A 的 C-alpha 接触图...")
    ca_coords, ca_residue_ids = get_coords(chain_a, "CA")
    if ca_coords.shape[0] > 0 : # 确保有坐标
        ca_contact_map = calculate_contact_map(ca_coords, CA_THRESHOLD, ca_residue_ids, pdb_id, "CA")
        # 你可以将 ca_contact_map (NumPy数组) 保存到文件，例如使用 np.save()
        if ca_contact_map is not None:
            np.save(os.path.join(OUTPUT_DIRECTORY, f"{pdb_id}_chainA_CA_contact_map.npy"), ca_contact_map)
    else:
        print(f"警告: 未能从 {pdb_id} Chain A 提取 C-alpha 坐标。")

    # --- C-beta 接触图 ---
    print(f"计算 {pdb_id} Chain A 的 C-beta 接触图...")
    cb_coords, cb_residue_ids = get_coords(chain_a, "CB") # get_coords 函数会处理Glycine
    if cb_coords.shape[0] > 0: # 确保有坐标
        cb_contact_map = calculate_contact_map(cb_coords, CB_THRESHOLD, cb_residue_ids, pdb_id, "CB")
        if cb_contact_map is not None:
            np.save(os.path.join(OUTPUT_DIRECTORY, f"{pdb_id}_chainA_CB_contact_map.npy"), cb_contact_map)
    else:
        print(f"警告: 未能从 {pdb_id} Chain A 提取 C-beta 坐标。")

    # --- 比较同一结构的 Calpha 和 Cbeta contact maps ---
    # 在这一步，你可以加载刚刚保存的 .npy 文件，或者直接使用上面得到的 ca_contact_map 和 cb_contact_map 变量
    # 进行比较。例如，计算它们之间的差异，或者并排显示它们。
    # 这里只作提示，具体比较方式取决于你的需求。
    if ca_contact_map is not None and cb_contact_map is not None:
        if ca_contact_map.shape == cb_contact_map.shape:
            diff_map = np.abs(ca_contact_map - cb_contact_map)
            common_contacts = np.logical_and(ca_contact_map, cb_contact_map).sum()
            ca_only_contacts = (ca_contact_map - np.logical_and(ca_contact_map, cb_contact_map)).sum()
            cb_only_contacts = (cb_contact_map - np.logical_and(ca_contact_map, cb_contact_map)).sum()
            print(f"{pdb_id}: 共同接触点数: {common_contacts}, "
                  f"仅Calpha接触点数: {ca_only_contacts}, "
                  f"仅Cbeta接触点数: {cb_only_contacts}")
        else:
            print(f"警告: {pdb_id} 的 CA 和 CB 接触图维度不匹配，无法直接比较。")
            print(f"CA shape: {ca_contact_map.shape}, CB shape: {cb_contact_map.shape}")

print("\n所有处理完成！")