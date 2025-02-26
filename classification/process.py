from Bio import SeqIO
from Bio.Seq import Seq
import collections
import csv
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
import os
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
os.chdir(script_dir)
# 配置：选择输出格式（'csv', 'tsv' 或 'xlsx'）
output_format = 'csv'  # 可改为 'tsv' 或 'xlsx'
output_file = 'sequence_matches'  # 输出文件基础名称（会自动添加扩展名）

# 加载基因组 FASTA 文件
fasta_file = "extracted.fa"  # 替换为你的 FASTA 文件路径
genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# 将基因组序列转换为大写，确保匹配时不区分大小写
for record in genome.values():
    record.seq = record.seq.upper()

# 从文本文件读取输入序列
input_file = "sequences.txt"  # 替换为你的输入序列文件路径
with open(input_file, "r") as f:
    queries = [line.strip().upper() for line in f if line.strip()]

# 定义函数：查找单个查询序列的匹配
def find_matches(query):
    query_seq = Seq(query)
    rev_comp_query = query_seq.reverse_complement()
    matches = []
    for chrom, record in genome.items():
        seq = record.seq
        # 在正链上查找所有匹配
        pos = -1
        while True:
            pos = seq.find(query_seq, pos + 1)
            if pos == -1:
                break
            start = pos + 1  # 转换为基于1的坐标
            end = pos + len(query)  # 基于1的结束坐标
            matches.append({
                "sequence": query,
                "chr": chrom,
                "start": start,
                "end": end,
                "strand": "+"
            })

        # 在负链（反向互补）上查找所有匹配
        pos = -1
        while True:
            pos = seq.find(rev_comp_query, pos + 1)
            if pos == -1:
                break
            start = pos + 1  # 转换为基于1的坐标
            end = pos + len(query)  # 基于1的结束坐标
            matches.append({
                "sequence": query,
                "chr": chrom,
                "start": start,
                "end": end,
                "strand": "-"
            })
    return matches

# 使用多进程并行搜索
if __name__ == '__main__':
    # 创建进程池，使用所有可用 CPU 核心
    with mp.Pool(mp.cpu_count()) as pool:
        # 并行执行 find_matches 函数，tqdm 显示进度
        results = list(tqdm(pool.imap(find_matches, queries), total=len(queries), desc="正在搜索序列"))

    # 将结果展平为单一列表
    matches = [match for sublist in results for match in sublist]

    # 计算每个序列的匹配次数
    sequence_counts = collections.Counter([match["sequence"] for match in matches])

    # 准备输出数据，添加 'unique' 字段
    headers = ["sequence", "chr", "start", "end", "strand", "unique"]
    data = []
    for match in matches:
        unique = "yes" if sequence_counts[match["sequence"]] == 1 else "no"
        data.append([
            match["sequence"],
            match["chr"],
            match["start"],
            match["end"],
            match["strand"],
            unique
        ])

    # 根据选择的格式写入文件
    if output_format == 'csv':
        with open(f"{output_file}.csv", "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(data)
        print(f"结果已写入 {output_file}.csv")

    elif output_format == 'tsv':
        with open(f"{output_file}.tsv", "w", newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(headers)
            writer.writerows(data)
        print(f"结果已写入 {output_file}.tsv")

    elif output_format == 'xlsx':
        df = pd.DataFrame(data, columns=headers)
        df.to_excel(f"{output_file}.xlsx", index=False)
        print(f"结果已写入 {output_file}.xlsx")

    else:
        raise ValueError("无效的 output_format。请选择 'csv'、'tsv' 或 'xlsx'。")