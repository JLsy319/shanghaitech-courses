import requests
import xml.etree.ElementTree as ET
import pandas as pd
import time

# NCBI E-utilities 基本 URL
EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# 你的邮箱地址 (NCBI 建议提供，以便在有问题时联系你)
EMAIL_ADDRESS = "liusy0319@foxmail.com" # 请替换成你自己的邮箱

# Coronaviridae 的 Taxon ID
TAXON_ID = "11118"

def search_coronaviruses_genomes():
    """
    使用 esearch 搜索所有基因组已测序的冠状病毒。
    """
    print(f"开始搜索 Taxon ID: {TAXON_ID} (Coronaviridae) 的完整基因组...")
    esearch_url = f"{EUTILS_BASE_URL}esearch.fcgi"
    params = {
        "db": "nuccore", # 核酸数据库
        "term": f"txid{TAXON_ID}[Organism:exp] AND (complete genome[Title] OR complete sequence[Title])", # 搜索词：冠状病毒科的完整基因组
        "retmode": "xml",
        "retmax": "100000", # 尝试获取大量记录，实际数量可能小于此值
        "email": EMAIL_ADDRESS,
        "tool": "MyCoronavirusProjectScript" # 自定义工具名称
    }

    try:
        response = requests.get(esearch_url, params=params)
        response.raise_for_status() # 如果请求失败则抛出异常
        root = ET.fromstring(response.content)
        id_list = [id_elem.text for id_elem in root.findall(".//Id")]
        print(f"找到 {len(id_list)} 个符合条件的核酸记录 ID。")
        return id_list
    except requests.exceptions.RequestException as e:
        print(f"Esearch 请求失败: {e}")
        return []
    except ET.ParseError as e:
        print(f"解析 Esearch XML 响应失败: {e}")
        print(f"响应内容: {response.text[:500]}...") # 打印部分响应内容以帮助调试
        return []

def fetch_virus_details(ids):
    """
    使用 efetch 获取指定 ID 列表的病毒名称和 Taxon ID。
    NCBI 的 efetch 对 GET 请求的 URL 长度有限制，所以我们需要分批处理大量的 ID。
    或者，我们可以使用 POST 请求。
    这里我们使用 esummary，它更适合获取摘要信息，包括 Taxon ID 和物种名。
    """
    print(f"开始获取 {len(ids)} 个记录的详细信息...")
    virus_data = []
    batch_size = 200 # 一次获取多少个记录的摘要信息

    for i in range(0, len(ids), batch_size):
        batch_ids = ids[i:i+batch_size]
        print(f"  正在处理批次 {i//batch_size + 1} (ID: {i+1} 到 {i+len(batch_ids)})...")
        esummary_url = f"{EUTILS_BASE_URL}esummary.fcgi"
        params = {
            "db": "nuccore",
            "id": ",".join(batch_ids),
            "retmode": "xml",
            "email": EMAIL_ADDRESS,
            "tool": "MyCoronavirusProjectScript"
        }

        try:
            response = requests.post(esummary_url, data=params) # 使用 POST 避免 URL 过长
            response.raise_for_status()
            root = ET.fromstring(response.content)

            for docsum in root.findall(".//DocSum"):
                try:
                    # Taxon ID 通常在 Item Name="TaxId" 中
                    taxid_elem = docsum.find("./Item[@Name='TaxId']")
                    taxid = taxid_elem.text if taxid_elem is not None else None

                    # 物种名称通常在 Item Name="Title" 或 Item Name="Organism" 中
                    # Title 字段可能包含更多信息，而 Organism 更直接
                    # 我们需要找到 <Item Name="Organism" Type="String">Organism Name</Item>
                    # 或者解析 <Item Name="Title" Type="String">Virus name complete genome</Item>
                    # 有时物种名在 <SourceDb>...</SourceDb><SourceId>...</SourceId><Item Name="Caption">ACCESSION</Item>
                    # <Item Name="Title">Virus full name strain etc, complete genome.</Item>
                    # <Item Name="Extra">gi|12345|gb|AY278741.1|</Item>
                    # <Item Name="Gi">12345</Item>
                    # <Item Name="CreateDate">2003/10/20</Item>
                    # <Item Name="UpdateDate">2018/09/10</Item>
                    # <Item Name="Flags">512</Item>
                    # <Item Name="TaxId">227859</Item>
                    # <Item Name="Length">29727</Item>
                    # <Item Name="Status">live</Item>
                    # <Item Name="ReplacedBy"/>
                    # <Item Name="Comment">...<Organism>SARS coronavirus</Organism>...</Item>
                    # <Item Name="AccessionVersion">AY278741.1</Item>
                    
                    # 从 Title 提取病毒名通常更可靠，因为它通常是标准格式
                    # 但我们需要的是物种名，而不是整个序列的标题
                    # 尝试从 <Item Name="Organism"> 获取
                    organism_elem = docsum.find("./Item[@Name='Organism']")
                    if organism_elem is not None:
                        name = organism_elem.text
                    else: # 如果直接的 Organism 字段没有，尝试从 Title 解析
                        title_elem = docsum.find("./Item[@Name='Title']")
                        name = title_elem.text.split(",")[0].split(" strain")[0].split(" isolate")[0].strip() if title_elem is not None else "Unknown"
                    
                    if taxid and name:
                        virus_data.append({"Name": name, "Taxon ID": taxid})
                    else:
                        print(f"  警告: 记录缺少 Name 或 Taxon ID。Title: {docsum.find('./Item[@Name=\'Title\']').text if docsum.find('./Item[@Name=\'Title\']') is not None else 'N/A'}")

                except Exception as e_parse:
                    print(f"  解析 DocSum 条目时出错: {e_parse}")
                    # 打印有问题的 DocSum 以便调试
                    # print(ET.tostring(docsum, encoding='unicode'))
            
            time.sleep(0.4) # 遵守 NCBI 的请求频率限制 (每秒不超过3次，如果没API Key)

        except requests.exceptions.RequestException as e:
            print(f"  Esummary/Efetch 请求失败 (批次 {i//batch_size + 1}): {e}")
            if response:
                print(f"  响应状态码: {response.status_code}")
                print(f"  响应内容: {response.text[:500]}...")
        except ET.ParseError as e:
            print(f"  解析 Esummary/Efetch XML 响应失败 (批次 {i//batch_size + 1}): {e}")
            if response:
                print(f"  响应内容: {response.text[:500]}...")
                
    # 去重，因为不同的基因组片段可能属于同一个病毒 Taxon ID
    # 我们期望的是每个唯一的 (Name, Taxon ID) 对
    if virus_data:
        df_temp = pd.DataFrame(virus_data)
        df_unique = df_temp.drop_duplicates(subset=["Name", "Taxon ID"]).reset_index(drop=True)
        print(f"获取到 {len(df_unique)} 个独特的病毒条目。")
        return df_unique.to_dict('records')
    else:
        print("未能获取到任何病毒数据。")
        return []

if __name__ == "__main__":
    # 1. 搜索基因组 ID
    genome_ids = search_coronaviruses_genomes()

    if genome_ids:
        # 2. 获取病毒详细信息
        virus_details_list = fetch_virus_details(genome_ids)

        if virus_details_list:
            # 3. 创建 DataFrame 并保存到 Excel
            df = pd.DataFrame(virus_details_list)
            
            # 确保列的顺序
            df = df[["Name", "Taxon ID"]]
            
            # 去除完全重复的行 (以防万一之前的去重逻辑不够完美)
            df.drop_duplicates(inplace=True)
            
            # 按病毒名称排序
            df.sort_values(by="Name", inplace=True)
            
            output_filename = "excel_file_1_coronaviruses.xlsx"
            try:
                df.to_excel(output_filename, index=False)
                print(f"成功将结果保存到 {output_filename}")
                print(f"总共找到 {len(df)} 种基因组已测序的冠状病毒。")
            except Exception as e:
                print(f"保存到 Excel 文件失败: {e}")
        else:
            print("未能获取病毒的详细信息，无法生成 Excel 文件。")
    else:
        print("未能找到任何基因组 ID，无法继续。")