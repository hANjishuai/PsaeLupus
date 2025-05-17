"""
B细胞表位处理器 - 命令行扩展版

新增功能
时间监控：关键步骤耗时统计
智能合并：可调节窗口的动态合并算法
质量控制：最终长度校验与过滤
功能特性
智能合并策略：邻近短片段自动合并，保持表位连续性
灵活参数控制：可调节最小长度、合并窗口等关键参数
高效IO处理：使用生成器处理大文件，内存安全
质量监控：输出统计信息，辅助结果验证

>sp|P12345|ANTIGEN_HUMAN  # Header格式不限
MGSSHHHHHHSQPM...         # 预测结果要求：
                          # - 表位区域用大写字母表示
                          # - 非表位区域小写
AAAVLKEPWAS[GHR]C         # 示例有效序列段落

Header                Epitope
sp|P12345|ANTIGEN_HUMAN  mgrahs
sp|P12345|ANTIGEN_HUMAN  epwasg
sp|Q67890|PROT_X         vklcst

输入输出说明
参数          格式要求                 示例值
--input-fasta  BepiPred3预测FASTA      predict_results.fasta
--min-length   最小表位长度(默认5)      5
--merge-window 合并窗口大小(默认5)      3
--output-tsv  输出TSV路径             epitopes.tsv

关键参数说明表
参数	类型	默认值	作用域	功能说明
--input-fasta	PATH	必填	process-epitopes	指定BepiPred3预测结果文件路径
--min-length	INT	5	process-epitopes	控制最终表位的最小氨基酸数
--merge-window	INT	5	process-epitopes	允许合并的间隔最大残基数
--output-tsv	TEXT	epitopes.tsv	process-epitopes	指定结果文件输出路径


执行示例
# 基本用法
python epitope_tool.py process-epitopes 
  --input-fasta predict_results.fasta 
  --output-tsv my_epitopes.tsv

# 扩展参数
python epitope_tool.py process-epitopes 
  --min-length 7          # 提高长度阈值
  --merge-window 3       # 缩小合并窗口
  --output-tsv strict_epitopes.tsv


"""
import os
import re
import time
import functools
import glob
import shlex  # 替代 click.utils.split_arg
from typing import List, Tuple
import click
from Bio import SeqIO

# ================= CLI命令组 =================
@click.group()
def cli():
    """B细胞表位处理工具集"""
    pass

# ================= 核心算法 =================
def time_logger(func):
    """执行时间追踪装饰器"""
    @functools.wraps(func)  # 保留原函数元数据
    def _wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        click.echo(f"⏱️ {func.__name__} 耗时: {elapsed:.2f}s")
        return result
    return _wrapper

def find_ep_regions(seq: str, min_len: int, window: int) -> List[Tuple[int, int]]:
    """智能表位检测与合并算法"""
    # 阶段1: 原始区域检测
    raw_regions = [(m.start(), m.end()) for m in re.finditer(r'[A-Z]+', seq)]
    
    # 阶段2: 动态窗口合并
    merged = []
    for curr_s, curr_e in raw_regions:
        if merged and (curr_s - merged[-1][1]) <= window:
            prev_s, prev_e = merged.pop()
            new_s = min(prev_s, curr_s)
            new_e = max(prev_e, curr_e)
            merged.append((new_s, new_e))
        else:
            merged.append((curr_s, curr_e))
    
    # 阶段3: 质量控制
    return [(s, e) for s, e in merged if (e - s) >= min_len]

# ================= 子命令实现 =================
@cli.command()
@click.option("--input-fasta", required=True, type=str,
             help="输入路径（支持通配符和带空格路径，用引号包裹）")
@click.option("--min-length", default=5, show_default=True)
@click.option("--merge-window", default=5, show_default=True)
@click.option("--output-tsv", default="epitopes.tsv")
@time_logger
def process_epitopes(input_fasta: str, min_length: int, 
                    merge_window: int, output_tsv: str):
    """
    处理BepiPred3预测结果
    - 单个文件：data/file.fasta
    - 通配符：data/*.fasta
    - 带空格路径："path/with spaces/*.fa"
    - 混合输入："file1.fa file2.fasta data/*.fasta"
    """
    # 步骤1: 安全拆分输入参数（兼容带空格路径）
    try:
        raw_patterns = shlex.split(input_fasta)  # 用shlex处理引号和空格
    except ValueError as e:
        click.secho(f"❌ 路径解析失败: {str(e)}", fg="red")
        raise click.Abort()

    # 步骤2: 扩展通配符并验证文件
    file_paths = []
    for pattern in raw_patterns:
        # 扩展通配符（支持递归**）
        expanded = glob.glob(pattern, recursive=True)
        
        if not expanded:
            click.secho(f"⚠️ 无匹配项: {pattern}", fg="yellow")
            continue
            
        for path in expanded:
            if not os.path.isfile(path):
                click.secho(f"⛔ 忽略非文件: {path}", fg="yellow")
                continue
            file_paths.append(os.path.abspath(path))  # 存储绝对路径

    # 最终检查
    if not file_paths:
        click.secho("❌ 错误：没有有效的输入文件", fg="red")
        raise click.Abort()

    # 显示处理文件列表
    click.secho("\n📂 将处理以下文件：", bold=True)
    for f in file_paths:
        click.echo(f"  - {f}")

    # 初始化统计
    total_ep = 0
    protein_count = 0
    
    with open(output_tsv, 'w') as fout:
        fout.write("Header\tEpitope\n")
        
        for path in file_paths:
            try:
                for record in SeqIO.parse(input_fasta, "fasta"):
                    protein_count += 1
                    seq = str(record.seq)
            
                    # 执行核心算法
                    regions = find_ep_regions(seq, min_length, merge_window)
            
                    # 写入结果
                    for s, e in regions:
                        epitope = seq[s:e].upper()
                        fout.write(f"{record.id}\t{epitope}\n")
                        total_ep += 1
            except Exception as e:
                click.secho(f"⚠️ 处理文件 {path} 时出错: {str(e)}", fg="red")
                continue
    
    # 输出统计报告
    click.echo(f"\n📊 处理报告:")
    click.echo(f"• 处理蛋白数: {protein_count}")
    click.echo(f"• 发现表位数: {total_ep}")
    click.echo(f"• 输出文件: {output_tsv}")

if __name__ == "__main__":
    cli()
