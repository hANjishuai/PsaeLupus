"""
PDB ID查询器 - 增强版

功能特性
并发查询：多线程加速UniProt到PDB的映射查询
智能重试：自动重试失败请求，提升鲁棒性
进度可视：实时显示处理进度与统计信息
灵活输入：支持大文件处理与多种ID格式

输入输出说明
参数          格式要求                 示例值
--input-tsv   输入TSV文件路径         uniprot_ids.tsv
--id-column   UniProt ID列名         uniprot_id
--output-tsv  输出TSV路径            pdb_mapping.tsv
--threads     并发线程数(默认4)       4
--delay       请求延迟(秒,默认0.2)    0.3

执行示例
# 基本用法
python pdb_fetcher.py fetch-pdb-ids 
  --input-tsv input_data.tsv 
  --id-column uniprot_accession 
  --output-tsv pdb_results.tsv

# 扩展参数
python pdb_fetcher.py fetch-pdb-ids 
  --threads 8                # 提高并发数
  --delay 0.5               # 增加请求间隔
  --output-tsv strict_results.tsv
"""
import csv
import time
import functools 
from typing import List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial, lru_cache 
import click
from tqdm import tqdm
from rcsbapi.search import search_attributes as attrs
import requests  # 新增导入

# ================= CLI命令组 =================
@click.group()
def cli():
    """PDB B细胞表位处理工具集"""
    pass

# ================= 配置装饰器 =================
def retry(max_retries=3, delay=1):
    """请求重试装饰器"""
    def decorator(func):
        @functools.wraps(func)  # 需要functools模块
        def wrapper(*args, **kwargs):
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    if attempt == max_retries - 1:
                        raise
                    time.sleep(delay * (attempt + 1))
            return None
        return wrapper
    return decorator

# ================= 核心算法 =================
@retry(max_retries=3, delay=1)
@lru_cache(maxsize=1024)
def fetch_pdb_ids_cache(uniprot_id: str) -> List[str]:
    """带缓存的PDB ID查询"""
    query = attrs.rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession == uniprot_id
    return list(query())

def process_batch(records: List[dict], id_col: str, pbar) -> List[Tuple[str, str]]:
    """批量处理函数"""
    results = []
    for row in records:
        uid = row[id_col].strip()
        if not uid:
            continue
        
        try:
            pdb_list = fetch_pdb_ids_cache(uid)
            results.extend([(uid, pdb) for pdb in pdb_list])
        except Exception as e:
            pbar.write(f"⚠️ 查询失败: {uid} ({str(e)})")
            results.append((uid, ""))
        
        pbar.update(1)
        time.sleep(float(pbar.delay))  # 动态延迟控制
    
    return results

# ================= 新增验证功能 =================
def filter_valid_prot_ids(input_file: str, output_file: str, max_workers: int = 8) -> None:
    """带并发验证的优化版本"""
    # 读取Prot_id列表
    with open(input_file, 'r') as f:
        prot_ids = list({line.strip() for line in f if line.strip()})

    valid_ids = []
    
    # 线程安全的验证函数
    def validate_id(prot_id: str) -> Tuple[str, bool]:
        url = f"https://files.rcsb.org/download/{prot_id}.pdb"
        try:
            response = requests.head(url, timeout=10, allow_redirects=True)
            return (prot_id, response.ok)
        except requests.exceptions.RequestException:
            return (prot_id, False)

    # 并发验证
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(validate_id, pid) for pid in prot_ids]
        
        # 使用tqdm进度条
        with tqdm(total=len(prot_ids), desc="🔄 验证进度", unit="ID") as pbar:
            for future in as_completed(futures):
                prot_id, is_valid = future.result()
                if is_valid:
                    valid_ids.append(prot_id)
                    pbar.set_postfix_str(f"✅ Valid: {prot_id}")
                else:
                    pbar.set_postfix_str(f"❌ Invalid: {prot_id}")
                pbar.update(1)

    # 保存结果
    with open(output_file, 'w') as f:
        f.write('\n'.join(sorted(valid_ids)))

    click.echo(f"\n📊 验证完成：原始ID数 {len(prot_ids)} → 有效ID数 {len(valid_ids)}")
    click.echo(f"💾 结果保存至：{output_file}")

# ================= CLI命令实现 =================
@cli.command()
@click.option("--input-tsv", required=True, type=click.Path(exists=True),
             help="输入TSV文件路径（须包含表头）")
@click.option("--id-column", required=True, 
             help="UniProt ID所在列名")
@click.option("--output-tsv", default="pdb_mapping.tsv",
             help="输出TSV文件路径")
@click.option("--threads", default=4, show_default=True,
             help="并发查询线程数")
@click.option("--delay", default=0.2, show_default=True,
             help="请求间隔时间（秒）")
def fetch_pdb_ids(input_tsv: str, id_column: str, 
                 output_tsv: str, threads: int, delay: float):
    """
    UniProt到PDB的批量映射查询
    
    输入文件要求：
    - TSV格式且包含表头
    - 指定列包含有效的UniProt ID
    - 支持ID重复（自动去重）
    """
    # 初始化统计
    stats = {
        'total': 0,
        'success': 0,
        'failed': 0,
        'pdb_count': 0
    }

    # 读取输入文件
    try:
        with open(input_tsv, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            unique_ids = list({row[id_column].strip() for row in reader if row[id_column].strip()})
            stats['total'] = len(unique_ids)
    except Exception as e:
        click.echo(f"❌ 文件读取失败: {str(e)}")
        return

    # 准备进度条
    with tqdm(
        total=stats['total'], 
        desc="🔄 处理进度", 
        unit="ID",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]"
    ) as pbar:
        pbar.delay = delay  # 附加延迟参数
        
        # 多线程处理
        results = []
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []
            batch_size = max(1, stats['total'] // (threads * 2))
            
            # 分批提交任务
            for i in range(0, stats['total'], batch_size):
                batch = [{'id': uid} for uid in unique_ids[i:i+batch_size]]
                future = executor.submit(
                    partial(process_batch, id_col='id', pbar=pbar),
                    batch
                )
                futures.append(future)
            
            # 收集结果
            for future in as_completed(futures):
                try:
                    batch_results = future.result()
                    results.extend(batch_results)
                except Exception as e:
                    pbar.write(f"❗ 批处理失败: {str(e)}")

        # 统计结果
        stats['success'] = len([r for r in results if r[1]])
        stats['failed'] = stats['total'] - stats['success']
        stats['pdb_count'] = len(results) - results.count(('', ''))

        # 写入输出文件
        try:
            with open(output_tsv, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(['UniProt_ID', 'PDB_ID'])
                writer.writerows(results)
        except Exception as e:
            click.echo(f"❌ 文件读取失败: {str(e)}")
            return

    # 输出统计报告
    click.echo("\n📊 执行报告:")
    click.echo(f"• 总处理ID数: {stats['total']}")
    click.echo(f"• 成功查询数: {stats['success']} ({stats['success']/stats['total']:.1%})")
    click.echo(f"• 失败查询数: {stats['failed']} ({stats['failed']/stats['total']:.1%})")
    click.echo(f"• 发现PDB总数: {stats['pdb_count']}")
    click.echo(f"• 输出文件: {output_tsv}")

# valid prot_ids
@cli.command()
@click.option("--input-file", required=True, type=click.Path(exists=True),
             help="待验证的Prot_id列表文件（每行一个ID）")
@click.option("--output-file", default="valid_prot_ids.txt",
             help="验证结果输出文件路径")
@click.option("--threads", default=8, show_default=True,
             help="并发验证线程数")
def validate_prot_ids(input_file: str, output_file: str, threads: int):
    """验证PDB ID有效性（通过HEAD请求检查文件存在性）"""
    try:
        filter_valid_prot_ids(
            input_file=input_file,
            output_file=output_file,
            max_workers=threads
        )
    except Exception as e:
        click.echo(f"❌ 验证过程出错: {str(e)}", err=True)
        
if __name__ == "__main__":
    cli()

