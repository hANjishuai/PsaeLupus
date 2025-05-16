# Get_Special_antigen_data.py
"""
# 下载FASTA（指定自定义查询和输出路径）
python skin_protein_tool.py fetch-fasta \
  --query "(reviewed:true) AND (organism_id:9606) AND (cc_tissue:skin)" \
  --output-fasta custom_skin_proteins.fasta

# 过滤元数据（指定输入文件）
python skin_protein_tool.py filter-data \
  --input-tsv uniprot_results.tsv \
  --filter-pattern "Secreted|Membrane" \
  --output-csv custom_filtered.csv

# 分割FASTA文件（输入文件需提前下载）
python skin_protein_tool.py split-fasta \
  --input-fasta antigen_proteins.fasta \
  --output-dir split_files \
  --chunk-size 50

设计特点
模块化设计：fetch_fasta 和 filter_data 作为独立子命令
参数灵活化：
可自定义查询语句、输入输出路径、过滤模式
默认值兼容原始需求 (skin_secreted_proteins.fasta 等)
进度可视化：click.echo 实时显示操作状态
错误处理：
HTTP 请求错误捕获
文件读取错误处理
编码安全：使用 urllib.parse.quote 处理特殊字符

输入输出说明
参数	格式要求	示例值
--query	UniProt 合法检索式	(reviewed:true) AND (organism_id:9606)
--input-tsv	UniProt 导出的 TSV 元数据文件	uniprot_results.tsv
--filter-pattern	Pandas 兼容的正则表达式	Secreted|Membrane|Nucleus
输出文件	FASTA/CSV 格式	.fasta, .csv

"""
import math 
import click
import requests
from pathlib import Path
import pandas as pd
from urllib.parse import quote

@click.group()
def cli():
    """UniProt 蛋白质数据获取工具"""
    pass

def get_fasta_uniport(query: str, output_fasta: str):
    """获取FASTA序列"""
    click.echo(f"⚡ 正在执行FASTA下载任务...")
    click.echo(f"📝 检索式: {query}")
    click.echo(f"💾 输出文件: {output_fasta}")
    
    encoded_query = quote(query)
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query={encoded_query}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(output_fasta, "w") as f:
            f.write(response.text)
        click.echo(f"✅ 成功下载 {len(response.text.split('>'))-1} 条FASTA序列")
    except Exception as e:
        click.echo(f"❌ 下载失败: {str(e)}", err=True)
        
def get_fasta_uniport2(query: str, output_fasta: str):
    """获取FASTA序列"""
    click.echo(f"⚡ 正在执行FASTA下载任务...")
    click.echo(f"📝 检索式: {query}")
    click.echo(f"💾 输出文件: {output_fasta}")
    url = query.replace("compressed=true", "compressed=false")
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        response.encoding = 'utf-8'
        with open(output_fasta, "w") as f:
            f.write(response.text)
        click.echo(f"✅ 成功下载 {len(response.text.split('>'))-1} 条FASTA序列")
    except Exception as e:
        click.echo(f"❌ 下载失败: {str(e)}", err=True)

def filter_metadata(input_tsv_url: str, filter_pattern: str, output_tsv: str, filted_output_tsv: str):
    """过滤元数据"""
    click.echo(f"\n⚡ 正在执行元数据过滤任务...")
    click.echo(f"📥 输入文件: {input_tsv_url}")
    click.echo(f"🔍 过滤模式: {filter_pattern}")
    click.echo(f"💾 输出文件: {output_tsv}")
    click.echo(f"💾 输出文件: {filted_output_tsv}")
    
    try:
        url = input_tsv_url.replace("compressed=true", "compressed=false")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        response.encoding = 'utf-8' 
        with open(output_tsv, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"✅ TSV 文件已保存至：{output_tsv}")

    except Exception as e:
        print(f"❌ 下载失败：{str(e)}")

    try:
        df = pd.read_csv(output_tsv, sep="\t")
        original_count = len(df)
        
        filtered_df = df[df["Subcellular location [CC]"].str.contains(filter_pattern, na=False)]
        filtered_count = len(filtered_df)
        
        filtered_df.to_csv(filted_output_tsv, index=False,sep="\t")
        click.echo(f"✅ 完成过滤: 原始数据 {original_count} 条 → 过滤后 {filtered_count} 条")
    except Exception as e:
        click.echo(f"❌ 过滤失败: {str(e)}", err=True)

def split_large_fasta(input_fasta: str, output_dir: str, chunk_size: int = 50):
    """修正版智能分割逻辑"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    base_name = Path(input_fasta).stem
    file_counter = 1
    current_chunk = []
    
    with open(input_fasta, 'r') as f:
        current_record = []
        for line in f:
            if line.startswith('>'):
                # 遇到新记录时处理当前缓存
                if current_record:
                    current_chunk.append(''.join(current_record))
                    current_record = []
                    
                    # 当积累到 chunk_size 时写入文件
                    if len(current_chunk) == chunk_size:
                        write_chunk(current_chunk, base_name, output_dir, file_counter)
                        file_counter += 1
                        current_chunk = []
                
                current_record.append(line)
            else:
                current_record.append(line)
        
        # 处理最后一个记录
        if current_record:
            current_chunk.append(''.join(current_record))
        
        # 写入剩余记录
        if current_chunk:
            write_chunk(current_chunk, base_name, output_dir, file_counter)

def write_chunk(records, base_name, output_dir, counter):
    """专用写入函数"""
    output_path = Path(output_dir) / f"{base_name}_part{counter}.fasta"
    with open(output_path, 'w') as out_f:
        out_f.write(''.join(records))
    click.echo(f"✂️ 生成分块文件: {output_path} (包含 {len(records)} 条序列)")


@cli.command()
@click.option("--query", required=True, help="UniProt检索式")
@click.option("--output-fasta", default="skin_secreted_proteins.fasta", 
             help="FASTA输出文件名")
def fetch_fasta(query, output_fasta):
    """下载FASTA序列"""
    get_fasta_uniport(query, output_fasta)

@cli.command()
@click.option("--input-tsv_url", required=True, help="元数据TSV文件路径")
@click.option("--filter-pattern", default="Secreted|Cell membrane", 
             help="亚细胞定位过滤正则")
@click.option("--output-tsv", default="proteins_meta.tsv",
             help="蛋白质序列meta_TSV文件名")
@click.option("--filted_output_tsv", default="filtered_proteins_meta.tsv",
             help="过滤后的蛋白质序列meta_TSV文件名")
def filter_data(input_tsv_url, filter_pattern, output_tsv, filted_output_tsv):
    """过滤蛋白质元数据"""
    filter_metadata(input_tsv_url, filter_pattern, output_tsv, filted_output_tsv)

@cli.command()
@click.option("--query", required=True, help="UniProt检索式")
@click.option("--output-fasta", default="skin_secreted_proteins.fasta", 
             help="FASTA输出文件名")
def fetch_fasta2(query, output_fasta):
    """下载FASTA序列"""
    get_fasta_uniport2(query, output_fasta)

@cli.command()
@click.option("--input-fasta", required=True, help="输入FASTA文件路径")
@click.option("--output-dir", default="split_results", help="分割文件输出目录")
@click.option("--chunk-size", default=50, help="每个文件最大序列数")
def split_fasta(input_fasta, output_dir, chunk_size):
    """优化后的分割命令"""
    click.echo(f"\n🔍 开始处理: {input_fasta}")
    
    if not Path(input_fasta).exists():
        click.echo(f"❌ 错误: 文件不存在", err=True)
        return
    
    with open(input_fasta, 'r') as f:
        seq_count = sum(1 for line in f if line.startswith('>'))
    
    click.echo(f"📊 总序列数: {seq_count}")
    
    if seq_count <= chunk_size:
        click.echo("✅ 无需分割")
        return
    
    click.echo(f"⚡ 开始分割（每 {chunk_size} 条/文件）...")
    split_large_fasta(input_fasta, output_dir, chunk_size)
    click.echo(f"🎉 分割完成！共生成 {math.ceil(seq_count/chunk_size)} 个文件")


if __name__ == "__main__":
    cli()

