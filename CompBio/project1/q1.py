from Bio import Entrez
import pandas as pd

Entrez.email = "justin_liu319@sina.com"

# 定义搜索函数
def fetch_coronaviridae_taxids():
    # 搜索冠状病毒科（Taxon ID: 11118）下所有子类
    handle = Entrez.esearch(db="taxonomy", term="txid11118[Subtree] AND genome sequenced[prop]", retmax=1000)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

# 获取所有Taxon ID
taxids = fetch_coronaviridae_taxids()

# 提取每个Taxon的详细信息
data = []
for taxid in taxids:
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    for record in records:
        name = record["ScientificName"]
        # 检查是否有基因组测序标记
        if "genbank" in record["OtherNames"]:
            data.append({"Name": name, "Taxon ID": taxid})
    handle.close()

# 保存到Excel
df = pd.DataFrame(data)
df.to_excel("project1/pj1-q1.xlsx", index=False)

