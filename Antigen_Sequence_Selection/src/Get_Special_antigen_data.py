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

import click
import requests
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

if __name__ == "__main__":
    cli()

