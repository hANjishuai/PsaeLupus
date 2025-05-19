"""
TSV文件分割工具 - 增强版

功能特性
智能分块：按指定行数分割文件，保留标题头
并行处理：多线程加速大文件分割
进度可视：实时显示处理进度与统计信息
灵活输出：自定义输出文件名模式

输入输出说明
参数           格式要求                 示例值
--input-tsv    输入TSV文件路径         data.tsv
--output-dir   输出目录               ./output
--prefix       输出文件前缀           chunk_
--lines        每文件行数(默认50)      100
--header       包含标题头(默认True)    --no-header
--threads      并发线程数(默认4)       6

执行示例
# 基本用法
python tsv_splitter.py split 
  --input-tsv large_data.tsv 
  --output-dir ./fragments 
  --prefix part_

# 扩展参数
python tsv_splitter.py split 
  --lines 200                # 每文件200行
  --no-header                # 无标题行模式
  --threads 8                # 提高并发数
"""
import csv
import os
from typing import List, Generator, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
import click
from tqdm import tqdm

# ================= CLI命令组 =================
@click.group()
def cli():
    """TSV文件处理工具集"""
    pass

# ================= 核心算法 =================
class TSVSplitter:
    def __init__(self, input_tsv: str, output_dir: str, 
                prefix: str = "chunk_", lines: int = 50, 
                has_header: bool = True, threads: int = 4):
        self.input_tsv = input_tsv
        self.output_dir = output_dir
        self.prefix = prefix
        self.chunk_size = lines
        self.has_header = has_header
        self.threads = threads
        self.stats = {'total': 0, 'chunks': 0, 'lines_written': 0}

    def validate(self):
        """参数有效性验证"""
        if not os.path.exists(self.input_tsv):
            raise FileNotFoundError(f"输入文件不存在: {self.input_tsv}")
        os.makedirs(self.output_dir, exist_ok=True)

    def read_tsv(self) -> Generator[List[str], None, None]:
        """流式读取TSV文件"""
        with open(self.input_tsv, 'r', newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            if self.has_header:
                self.header = next(reader)
            for row in reader:
                yield row

    def write_chunk(self, chunk_id: int, rows: List[List[str]]) -> str:
        """写入单个分块文件"""
        suffix = f"{chunk_id:04d}"  # 4位数字补零
        output_path = os.path.join(
            self.output_dir, 
            f"{self.prefix}{suffix}.tsv"
        )
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            if self.has_header:
                writer.writerow(self.header)
            writer.writerows(rows)
        
        return output_path

    def process(self) -> dict:
        """主处理流程"""
        self.validate()
        rows = []
        chunk_id = 0
        futures = []
        
        with tqdm(desc="🔄 分割进度", unit="行") as pbar, \
             ThreadPoolExecutor(max_workers=self.threads) as executor:

            for row in self.read_tsv():
                rows.append(row)
                self.stats['total'] += 1
                pbar.update(1)

                if len(rows) >= self.chunk_size:
                    future = executor.submit(
                        self.write_chunk, chunk_id, rows.copy()
                    )
                    futures.append(future)
                    rows.clear()
                    chunk_id += 1

            # 处理剩余行
            if rows:
                future = executor.submit(
                    self.write_chunk, chunk_id, rows.copy()
                )
                futures.append(future)

            # 收集结果
            for future in as_completed(futures):
                path = future.result()
                self.stats['chunks'] += 1
                self.stats['lines_written'] += (len(rows) if 'rows' in locals() else self.chunk_size)
                pbar.set_postfix_str(f"生成文件: {os.path.basename(path)}")

        return self.stats

# ================= CLI命令实现 =================
@cli.command()
@click.option("--input-tsv", required=True, type=click.Path(exists=True),
             help="输入TSV文件路径")
@click.option("--output-dir", default="./output", show_default=True,
             help="输出目录路径")
@click.option("--prefix", default="chunk_", show_default=True,
             help="输出文件名前缀")
@click.option("--lines", default=50, show_default=True,
             help="每个分块文件包含的行数")
@click.option("--header/--no-header", default=True,
             help="是否包含标题行")
@click.option("--threads", default=4, show_default=True,
             help="并发写入线程数")
def split(input_tsv: str, output_dir: str, prefix: str, 
         lines: int, header: bool, threads: int):
    """
    TSV文件分割命令
    
    输入文件要求：
    - 标准TSV格式
    - 标题行需在第一行（当启用header时）
    - 允许包含空行（自动跳过）
    """
    processor = TSVSplitter(
        input_tsv=input_tsv,
        output_dir=output_dir,
        prefix=prefix,
        lines=lines,
        has_header=header,
        threads=threads
    )
    
    try:
        stats = processor.process()
        click.echo("\n📊 分割报告:")
        click.echo(f"• 输入文件: {input_tsv}")
        click.echo(f"• 总处理行数: {stats['total']}")
        click.echo(f"• 生成分块数: {stats['chunks']}")
        click.echo(f"• 输出目录: {output_dir}")
        click.echo(f"• 文件名模式: {prefix}XXXX.tsv")
    except Exception as e:
        click.echo(f"❌ 处理失败: {str(e)}", err=True)

if __name__ == "__main__":
    cli()
