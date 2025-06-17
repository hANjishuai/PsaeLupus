"""
抗原分析工具集 - 整合版

子命令列表
calculate-ratios : 分块计算表位比率

"""

import sys
import re
import csv
import time
import click
import functools
from pathlib import Path
import matplotlib.pyplot as plt
from typing import Dict, Generator, Tuple , DefaultDict, Set , List , Optional
from collections import defaultdict
import pandas as pd
import plotly.express as px
import plotly.io as pio


# ================= 基础架构 =================
class BaseCommand:
    """命令基类"""
    def __init__(self, input_path: Path, output_path: Path):
        self.input = input_path
        self.output = output_path
        
    def validate(self):
        """输入验证"""
        if not self.input.exists():
            raise FileNotFoundError(f"输入路径不存在: {self.input}")
        self.output.parent.mkdir(parents=True, exist_ok=True)

# ================= 比率计算模块 =================
class RatioCalculator(BaseCommand):
    """分块计算表位比率"""
    def __init__(self, input_dir: Path, output_file: Path, 
                batch_size: int = 50, min_residues: int = 50):
        super().__init__(input_dir, output_file)
        self.batch_size = batch_size
        self.min_residues = min_residues
        self.timer = Timer()

    def process_single_file(self, file_path: Path) -> Dict[str, Tuple[int, int]]:
        """处理单个文件"""
        counts = defaultdict(lambda: {'total': 0, 'epitope': 0})
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f)
                if not {'pdb', 'epitope'}.issubset(reader.fieldnames):
                    raise ValueError("缺少必需列")
                    
                for row in reader:
                    antigen = row['pdb']
                    is_epitope = row['epitope'].lower() in ['true', '1', 'yes']
                    
                    counts[antigen]['total'] += 1
                    if is_epitope:
                        counts[antigen]['epitope'] += 1
        except Exception as e:
            click.secho(f"❌ 处理失败: {file_path} - {str(e)}", fg='red')
        return counts

    def batch_processor(self, all_files: list) -> Generator[Tuple[str, Dict], None, None]:
        """分块处理文件"""
        self.timer.start('processing')
        for i, file_path in enumerate(all_files, 1):
            counts = self.process_single_file(file_path)
            yield from counts.items()
            
            # 进度显示
            if i % 10 == 0:
                click.secho(f"🔄 已处理 {i}/{len(all_files)} 文件", fg='green')
        self.timer.end('processing')

    def execute(self):
        """执行入口"""
        # 获取文件列表
        all_files = list(self.input.glob("*.csv"))
        if not all_files:
            raise FileNotFoundError("目录中没有CSV文件")
        
        click.secho(f"📂 发现 {len(all_files)} 个抗原文件", fg='blue')
        
        # 初始化结果存储
        global_counts = defaultdict(lambda: {'total': 0, 'epitope': 0})
        
        # 分块处理
        with open(self.output, 'w', newline='') as f_out:
            writer = csv.DictWriter(f_out, fieldnames=[
                'Antigen', 'Total', 'Epitope', 'Ratio'
            ])
            writer.writeheader()
            
            for batch_start in range(0, len(all_files), self.batch_size):
                batch_files = all_files[batch_start:batch_start+self.batch_size]
                
                # 处理当前批次
                for antigen, counts in self.batch_processor(batch_files):
                    global_counts[antigen]['total'] += counts['total']
                    global_counts[antigen]['epitope'] += counts['epitope']
                
                # 写入并清空缓存
                for antigen, data in global_counts.items():
                    if data['total'] < self.min_residues:
                        continue
                    
                    ratio = data['epitope'] / data['total'] if data['total'] > 0 else 0.0
                    writer.writerow({
                        'Antigen': antigen,
                        'Total': data['total'],
                        'Epitope': data['epitope'],
                        'Ratio': f"{ratio:.4f}"
                    })
                global_counts.clear()

# ================= 比率直方图生成器 =================
class RatioHistogramPlotter(BaseCommand):
    """比率直方图生成器"""
    def __init__(self, input_file: Path, output_image: Path, 
                bins: int = 20, color: str = 'skyblue'):
        super().__init__(input_file, output_image)
        self.bins = bins
        self.color = color
        self.timer = Timer()

    def load_ratios(self) -> list:
        """加载比率数据"""
        ratios = []
        with open(self.input, 'r') as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                try:
                    ratio = float(row['Ratio'])
                    ratios.append(ratio)
                except (KeyError, ValueError) as e:
                    click.secho(f"⚠️ 无效数据行: {row} - {str(e)}", fg='yellow')
        return ratios

    def execute(self):
        """生成直方图"""
        # 加载数据
        self.timer.start('data_loading')
        click.secho("📊 正在加载数据...", fg='blue')
        ratios = self.load_ratios()
        if not ratios:
            raise ValueError("未找到有效比率数据")
        self.timer.end('data_loading')

        # 创建图表
        self.timer.start('plotting')
        plt.figure(figsize=(10, 6))
        n, bins, patches = plt.hist(
            ratios, 
            bins=self.bins,
            color=self.color,
            edgecolor='black',
            alpha=0.7
        )

        # 样式设置
        plt.title('Histogram of Antigen Epitope Ratio Distribution', fontsize=14)
        plt.xlabel('Epitope Proportion %', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.grid(axis='y', alpha=0.75)
        
        # 保存输出
        plt.tight_layout()
        plt.savefig(self.output, dpi=300)
        self.timer.end('plotting')

        click.secho(f"✅ 图表已保存至 {self.output}", fg='green')

# ================= 比率直方图+拟合曲线生成器 =================
class EnhancedRatioPlotter(RatioHistogramPlotter):
    """支持拟合曲线的直方图"""
    def __init__(self, input_file: Path, output_image: Path, 
                bins: int = 20, color: str = 'skyblue',
                fit_method: str = 'kde', 
                curve_color: str = 'darkred',
                show_legend: bool = True):
        super().__init__(input_file, output_image, bins, color)
        self.fit_method = fit_method
        self.curve_color = curve_color
        self.show_legend = show_legend
        self._validate_fit_method()

    def _validate_fit_method(self):
        """验证拟合方法有效性"""
        valid_methods = ['kde', 'normal', 'none']
        if self.fit_method not in valid_methods:
            raise ValueError(f"无效拟合方法，可选：{valid_methods}")

    def _add_fit_curve(self, ratios: list):
        """添加拟合曲线"""
        import numpy as np
        from scipy.stats import gaussian_kde, norm

        # 生成统一x轴范围
        xmin, xmax = min(ratios), max(ratios)
        x = np.linspace(xmin, xmax, 500)

        # 选择拟合方法
        if self.fit_method == 'kde':
            try:
                kde = gaussian_kde(ratios)
                y = kde(x)
                label = 'Kernel Density Estimation(KDE)'
            except Exception as e:
                click.secho(f"⚠️ KDE拟合失败：{str(e)} 改用正态分布", fg='yellow')
                self.fit_method = 'normal'

        if self.fit_method == 'normal':
            mu, std = np.mean(ratios), np.std(ratios)
            y = norm.pdf(x, mu, std)
            label = f'Normal Distribution (μ={mu:.2f}, σ={std:.2f})'

        if self.fit_method != 'none':
            plt.plot(x, y, color=self.curve_color, lw=2, label=label)
            if self.show_legend:
                plt.legend()

    def execute(self):
        """生成增强版图表"""
        # 加载基础数据
        ratios = self.load_ratios()
        if not ratios:
            raise ValueError("无有效数据用于拟合")

        # 创建基础直方图
        plt.figure(figsize=(10, 6))
        n, bins, patches = plt.hist(
            ratios, 
            bins=self.bins,
            color=self.color,
            edgecolor='black',
            alpha=0.7,
            density=True  # 标准化为概率密度
        )

        # 添加拟合曲线
        if self.fit_method != 'none':
            self.timer.start('curve_fitting')
            self._add_fit_curve(ratios)
            self.timer.end('curve_fitting')

        # 增强样式
        plt.title('Distribution of Antigen Ratios and Fitting Curves', fontsize=14)
        plt.xlabel('Epitope Proportion', fontsize=12)
        plt.ylabel('Probability Density', fontsize=12)
        plt.grid(axis='y', alpha=0.75)

        # 保存输出
        plt.tight_layout()
        plt.savefig(self.output, dpi=300)
        plt.close()

# ================= 增强数据整合模块 =================
class EnhancedDataMerger:
    """整合PDB映射、蛋白元数据和表位比率（保留所有表位数据）"""
    def __init__(self, pdb_mapping: Path, proteins_meta: Path, epitope_ratios: Path, output: Path):
        self.pdb_mapping = pdb_mapping
        self.proteins_meta = proteins_meta
        self.epitope_ratios = epitope_ratios
        self.output = output
        self.timer = Timer()
        
    def validate(self):
        """验证所有输入文件存在"""
        required_files = [
            (self.pdb_mapping, "PDB映射文件"),
            (self.proteins_meta, "蛋白元数据文件"),
            (self.epitope_ratios, "表位比率文件")
        ]
        
        for path, name in required_files:
            if not path.exists():
                raise FileNotFoundError(f"{name}不存在: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)
    
    def load_pdb_mapping(self) -> Dict[str, Set[str]]:
        """加载PDB到UniProt的映射"""
        mapping = defaultdict(set)
        with open(self.pdb_mapping, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                pdb_id = row['PDB_ID'].strip()
                uniprot_id = row['UniProt_ID'].strip()
                if pdb_id and uniprot_id:
                    mapping[pdb_id].add(uniprot_id)
        return mapping
    
    def load_proteins_meta(self) -> Dict[str, Dict]:
        """加载蛋白元数据，建立UniProt ID到元数据的映射"""
        meta = {}
        with open(self.proteins_meta, 'r', encoding='utf-8') as f:
            # 处理可能的列名变化
            fieldnames = next(csv.reader(f, delimiter='\t'))
            f.seek(0)
            
            # 查找Tissue specificity列
            tissue_col = None
            for col in fieldnames:
                if 'Tissue specificity' in col:
                    tissue_col = col
                    break
            
            reader = csv.DictReader(f, delimiter='\t', fieldnames=fieldnames)
            next(reader)  # 跳过标题行
            
            for row in reader:
                uniprot_id = row['Entry'].strip()
                if not uniprot_id:
                    continue
                
                entry_name = row.get('Entry Name', '').strip()
                tissue_specificity = row.get(tissue_col, '').strip() if tissue_col else ''
                
                # 清理PDB字段（可能包含多个PDB ID）
                pdb_ids = []
                if 'PDB' in row and row['PDB']:
                    pdb_list = row['PDB'].split(';')
                    for pdb in pdb_list:
                        pdb = pdb.strip()
                        if pdb:
                            pdb_ids.append(pdb.split()[0])  # 取PDB ID部分，忽略可能的注释
                
                meta[uniprot_id] = {
                    'Entry_Name': entry_name,
                    'Tissue_specificity': tissue_specificity,
                    'PDB_IDs': set(pdb_ids)
                }
        return meta
    
    def split_entry_name(self, entry_name: str) -> Tuple[str, str]:
        """拆分Entry Name为Gene和Species"""
        if entry_name == "NA" or not entry_name:
            return "NA", "NA"
        
        # 查找最后一个下划线的位置
        last_underscore = entry_name.rfind('_')
        
        if last_underscore == -1:
            # 没有下划线，整个作为gene
            return entry_name, "UNKNOWN"
        
        gene = entry_name[:last_underscore]
        species = entry_name[last_underscore + 1:]
        
        # 清理物种名称（如HUMAN -> Human）
        species_clean = species.capitalize()
        if species_clean == "Human":
            species_clean = "Homo sapiens"
        elif species_clean == "Mouse":
            species_clean = "Mus musculus"
        elif species_clean == "Rat":
            species_clean = "Rattus norvegicus"
        
        return gene, species_clean
    
    def merge_data(self):
        """执行数据整合 - 保留所有epitope_ratios数据"""
        self.timer.start('overall')
        
        # 步骤1: 加载数据
        self.timer.start('loading_data')
        pdb_to_uniprot = self.load_pdb_mapping()
        uniprot_to_meta = self.load_proteins_meta()
        self.timer.end('loading_data')
        click.secho(f"📊 加载 {len(pdb_to_uniprot)} 个PDB映射和 {len(uniprot_to_meta)} 个蛋白元数据", fg='blue')
        
        # 步骤2: 处理表位比率数据
        self.timer.start('merging')
        results = []
        missing_pdb = set()
        missing_uniprot = set()
        found_count = 0
        
        with open(self.epitope_ratios, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            total_rows = sum(1 for _ in reader)  # 获取总行数用于进度显示
            f.seek(0)
            next(reader)  # 跳过标题行
            
            for i, row in enumerate(reader, 1):
                antigen = row['Antigen']
                # 进度显示
                if i % 100 == 0 or i == total_rows:
                    click.secho(f"🔄 处理中: {i}/{total_rows} 行 ({i/total_rows:.1%})", fg='cyan')
                
                # 解析抗原信息
                parts = antigen.split('_')
                pdb_id = parts[0] if parts else ""
                chain_id = parts[1] if len(parts) > 1 else ""
                
                # 准备结果行（保留所有原始数据）
                result_row = {
                    'PDB_ID': pdb_id,
                    'Chain_ID': chain_id,
                    'Uniprot_ID': "NA",
                    'Entry_Name': "NA",
                    'Gene': "NA",
                    'Species': "NA",
                    'Tissue_specificity': "NA",
                    'Total': row['Total'],
                    'Epitope': row['Epitope'],
                    'Ratio': row['Ratio'],
                    'Ratio_float': 0.0  # 用于排序
                }
                
                # 尝试通过PDB映射找到UniProt ID
                uniprot_ids = pdb_to_uniprot.get(pdb_id, set())
                found = False
                
                for uniprot_id in uniprot_ids:
                    if uniprot_id in uniprot_to_meta:
                        meta = uniprot_to_meta[uniprot_id]
                        result_row['Uniprot_ID'] = uniprot_id
                        result_row['Entry_Name'] = meta['Entry_Name']
                        result_row['Tissue_specificity'] = meta['Tissue_specificity']
                        
                        # 拆分Entry Name
                        gene, species = self.split_entry_name(meta['Entry_Name'])
                        result_row['Gene'] = gene
                        result_row['Species'] = species
                        
                        found = True
                        found_count += 1
                        break  # 使用第一个匹配项
                
                # 如果通过PDB映射没找到，尝试通过PDB列表查找
                if not found:
                    for uniprot_id, meta in uniprot_to_meta.items():
                        if pdb_id in meta['PDB_IDs']:
                            result_row['Uniprot_ID'] = uniprot_id
                            result_row['Entry_Name'] = meta['Entry_Name']
                            result_row['Tissue_specificity'] = meta['Tissue_specificity']
                            
                            # 拆分Entry Name
                            gene, species = self.split_entry_name(meta['Entry_Name'])
                            result_row['Gene'] = gene
                            result_row['Species'] = species
                            
                            found = True
                            found_count += 1
                            break
                
                # 记录缺失信息
                if not found:
                    if uniprot_ids:
                        missing_uniprot.update(uniprot_ids)
                    elif pdb_id:
                        missing_pdb.add(pdb_id)
                
                # 转换Ratio为浮点数用于排序
                try:
                    result_row['Ratio_float'] = float(row['Ratio'])
                except ValueError:
                    result_row['Ratio_float'] = 0.0
                
                results.append(result_row)
        
        self.timer.end('merging')
        
        # 步骤3: 按Ratio降序排序
        self.timer.start('sorting')
        click.secho("🔢 按Ratio降序排序...", fg='blue')
        results.sort(key=lambda x: x['Ratio_float'], reverse=True)
        self.timer.end('sorting')
        
        # 步骤4: 写入结果
        self.timer.start('writing')
        with open(self.output, 'w', newline='', encoding='utf-8') as f_out:
            fieldnames = [
                'PDB_ID', 'Chain_ID', 'Uniprot_ID', 'Gene', 'Species', 
                'Tissue_specificity', 'Total', 'Epitope', 'Ratio'
            ]
            writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            for row in results:
                # 准备最终输出行（移除临时字段）
                output_row = {k: v for k, v in row.items() if k != 'Ratio_float' and k != 'Entry_Name'}
                writer.writerow(output_row)
        
        # 输出统计信息
        click.secho(f"✅ 成功整合 {len(results)} 条记录", fg='green')
        click.secho(f"🔍 匹配到元数据的记录: {found_count} ({found_count/len(results):.1%})", fg='green')
        click.secho(f"📊 最高Ratio: {results[0]['Ratio']} (PDB: {results[0]['PDB_ID']})", fg='green')
        click.secho(f"📊 最低Ratio: {results[-1]['Ratio']} (PDB: {results[-1]['PDB_ID']})", fg='green')
        
        if missing_pdb:
            click.secho(f"⚠️ 警告: {len(missing_pdb)} 个PDB ID未找到映射: {', '.join(sorted(missing_pdb)[:5])}{'...' if len(missing_pdb) > 5 else ''}", fg='yellow')
        if missing_uniprot:
            click.secho(f"⚠️ 警告: {len(missing_uniprot)} 个UniProt ID未找到元数据: {', '.join(sorted(missing_uniprot)[:5])}{'...' if len(missing_uniprot) > 5 else ''}", fg='yellow')
        
        self.timer.end('writing')
        self.timer.end('overall')
        
# ================= 增强版TreeMap可视化工具 =================
class TreeMapVisualizer:
    """TreeMap可视化工具 - 确保正确渲染矩形"""
    def __init__(self, input_file: Path, output_file: Path, 
                species_filter: Optional[str] = None, 
                min_ratio: float = 0.0,
                output_format: str = "png"):
        self.input_file = input_file
        self.output_file = output_file
        self.species_filter = species_filter
        self.min_ratio = min_ratio
        self.output_format = output_format.lower()
    
    def validate(self):
        """验证输入文件"""
        if not self.input_file.exists():
            raise FileNotFoundError(f"输入文件不存在: {self.input_file}")
        
        valid_formats = ["png", "svg", "jpg", "jpeg", "pdf"]
        if self.output_format not in valid_formats:
            raise ValueError(f"无效输出格式: {self.output_format}. 支持的格式: {', '.join(valid_formats)}")
        
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
    
    def load_and_filter_data(self) -> pd.DataFrame:
        """加载并过滤数据"""
        click.secho(f"📂 加载数据: {self.input_file}", fg='blue')
        df = pd.read_csv(self.input_file, sep='\t')
        
        # 基本数据验证
        required_columns = ['PDB_ID', 'Chain_ID', 'Gene', 'Species', 'Ratio']
        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            raise ValueError(f"输入文件缺少必要列: {', '.join(missing)}")
        
        # 转换比率列为数值类型
        df['Ratio'] = pd.to_numeric(df['Ratio'], errors='coerce')
        df.dropna(subset=['Ratio'], inplace=True)
        
        # 应用过滤器
        initial_count = len(df)
        
        if self.species_filter:
            click.secho(f"🔍 应用物种过滤: {self.species_filter}", fg='cyan')
            df = df[df['Species'].str.contains(self.species_filter, case=False, na=False)]
        
        if self.min_ratio > 0:
            click.secho(f"🔍 应用最小比率过滤: ≥{self.min_ratio:.2f}", fg='cyan')
            df = df[df['Ratio'] >= self.min_ratio]
        
        filtered_count = len(df)
        click.secho(f"✅ 成功加载 {initial_count} 条记录", fg='green')
        click.secho(f"🔍 过滤后保留 {filtered_count} 条记录 ({filtered_count/initial_count:.1%})", fg='cyan')
        
        return df
    
    def create_treemap(self, df: pd.DataFrame):
        """创建TreeMap可视化 - 简化方法确保正确渲染"""
        # 创建复合ID
        df['Composite_ID'] = df['PDB_ID'] + '_' + df['Chain_ID']
        
        # 创建层级数据
        df['Gene_PDB'] = df['Gene'] + '|' + df['Composite_ID']
        
        # 确保Ratio值大于0
        df = df[df['Ratio'] > 0]
        
        # 创建TreeMap - 使用plotly的简单方法
        fig = px.treemap(
            df,
            path=['Gene', 'Composite_ID'],
            values='Ratio',
            color='Gene',
            color_discrete_sequence=px.colors.qualitative.Dark24,
            title='B细胞表位预测TreeMap' + 
                "\n过滤阈值：" + str(self.min_ratio) + 
                '\n物种：' + self.species_filter,
            branchvalues='total',
            hover_name='Composite_ID',
            hover_data={
                'PDB_ID': True,
                'Chain_ID': True,
                'Ratio': ':.2%',
                'Species': True
            },
            height=1200,
            width=1600
        )
        
        # 更新布局确保正确渲染
        fig.update_layout(
            margin=dict(t=100, l=10, r=10, b=10),
            font=dict(
                family="Arial, sans-serif",
                size=14,
                color="black"
            ),
            paper_bgcolor='white',
            plot_bgcolor='white',
            uniformtext=dict(
                minsize=12,
                mode='hide'
            )
        )
        
        # 更新轨迹设置
        fig.update_traces(
            textinfo="label+value",
            texttemplate='<b>%{label}</b><br>%{value:.2%}',
            textfont=dict(
                size=14,
                color='white'
            ),
            textposition="middle center",
            marker_line_color='rgba(0,0,0,0.8)',
            marker_line_width=1.5,
            tiling=dict(
                packing='squarify',
                pad=10,
                squarifyratio=1
            )
        )
        
        return fig
    
    def save_output(self, fig):
        """保存输出文件 - 优化静态导出设置"""
        # 确保文件扩展名正确
        output_path = self.output_file.with_suffix(f".{self.output_format}")
        
        try:
            # 静态图片导出
            if pio.kaleido.scope is None:
                click.secho("❌ 导出静态图片需要安装kaleido包", fg='red')
                click.secho("💡 请运行: pip install kaleido", fg='yellow')
                sys.exit(1)
            
            # 设置高分辨率导出
            pio.kaleido.scope.default_format = self.output_format
            pio.kaleido.scope.default_width = 2000
            pio.kaleido.scope.default_height = 1500
            pio.kaleido.scope.default_scale = 3
            
            # 保存图片
            fig.write_image(output_path, engine="kaleido")
            click.secho(f"🖼️ {self.output_format.upper()}图片已保存至: {output_path}", fg='green')
        
        except Exception as e:
            click.secho(f"❌ 保存文件失败: {str(e)}", fg='red')
            if "kaleido" in str(e).lower():
                click.secho("💡 请确保已安装kaleido: pip install kaleido", fg='yellow')
            else:
                import traceback
                click.secho(traceback.format_exc(), fg='red')
   
    def execute(self):
        """执行可视化流程"""
        self.validate()
        df = self.load_and_filter_data()    
        if df.empty:
            click.secho("⚠️ 过滤后无数据可显示，请调整过滤条件", fg='yellow')
            return
        
        # 确保有足够的数据点
        if len(df) < 5:
            click.secho("⚠️ 数据点过少，无法生成有效的TreeMap", fg='yellow')
            return
        
        # 创建并保存TreeMap
        fig = self.create_treemap(df)
        self.save_output(fig)

# ================= CLI命令注册 =================
class Timer:
    """计时工具"""
    def __init__(self):
        self.stages = {}
    
    def start(self, name):
        self.stages[name] = time.perf_counter()
    
    def end(self, name):
        elapsed = time.perf_counter() - self.stages[name]
        click.secho(f"⏱️ {name} 耗时: {elapsed:.2f}s", fg='cyan')

def common_params(func):
    """共享参数"""
    @click.option("--input-dir", "-i", type=Path, required=True, 
                 help="输入目录路径")
    @click.option("--output", "-o", type=Path, default="ratios.tsv",
                 help="输出文件路径")
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper

@click.group()
def cli():
    """抗原数据分析工具"""
    pass

# ================= 比率计算模块 =================
@cli.command()
@common_params
@click.option("--batch-size", type=int, default=50,
             help="分块处理大小")
@click.option("--min-residues", type=int, default=50,
             help="最小残基数阈值")
def calculate_ratios(input_dir: Path, output: Path, batch_size: int, min_residues: int):
    """分块计算表位比率"""
    processor = RatioCalculator(
        input_dir=input_dir,
        output_file=output,
        batch_size=batch_size,
        min_residues=min_residues
    )
    
    try:
        processor.validate()
        click.secho("🚀 开始处理...", fg='green')
        processor.execute()
        click.secho(f"✅ 完成！结果保存至 {output}", fg='green', bold=True)
    except Exception as e:
        click.secho(f"❌ 错误: {str(e)}", fg='red')

#  ================= 比率直方图生成器 =================
@cli.command()
@click.option("--input", "-i", type=Path, required=True,
             help="比率计算结果文件路径（TSV格式）")
@click.option("--output", "-o", type=Path, default="ratio_histogram.png",
             help="输出图片路径（支持PNG/PDF/SVG）")
@click.option("--bins", type=int, default=20,
             help="直方图分箱数量")
@click.option("--color", default="skyblue",
             help="柱状图填充颜色")
def plot_histogram(input: Path, output: Path, bins: int, color: str):
    """生成抗原比率分布直方图"""
    plotter = RatioHistogramPlotter(
        input_file=input,
        output_image=output,
        bins=bins,
        color=color
    )
    
    try:
        plotter.validate()
        click.secho("🎨 正在生成直方图...", fg='blue')
        plotter.execute()
    except Exception as e:
        click.secho(f"❌ 绘图失败: {str(e)}", fg='red')

# ================= 比率直方图+拟合曲线生成器 =================
@cli.command()
@click.option("--input", "-i", type=Path, required=True,
             help="比率计算结果文件路径")
@click.option("--output", "-o", type=Path, default="enhanced_plot.png",
             help="输出图片路径")
@click.option("--bins", type=int, default=20,
             help="直方图分箱数")
@click.option("--color", default="skyblue",
             help="直方图颜色")
@click.option("--fit-method", type=click.Choice(['kde', 'normal', 'none']),
             default='kde', help="拟合方法：kde(默认)/normal/none")
@click.option("--curve-color", default="darkred",
             help="拟合曲线颜色")
@click.option("--show-legend/--no-legend", default=True,
             help="是否显示图例")
def enhanced_plot(input: Path, output: Path, bins: int, color: str,
                 fit_method: str, curve_color: str, show_legend: bool):
    """生成带拟合曲线的增强版图表"""
    plotter = EnhancedRatioPlotter(
        input_file=input,
        output_image=output,
        bins=bins,
        color=color,
        fit_method=fit_method,
        curve_color=curve_color,
        show_legend=show_legend
    )
    
    try:
        plotter.validate()
        click.secho("📈 正在生成增强图表...", fg='blue')
        plotter.execute()
        click.secho(f"✅ 图表保存至：{output}", fg='green')
    except Exception as e:
        click.secho(f"❌ 生成失败：{str(e)}", fg='red')

# ================= 增强数据整合命令 =================
@cli.command()
@click.option("--pdb-mapping", type=Path, required=True,
             help="PDB到UniProt的映射文件路径 (TSV格式)")
@click.option("--proteins-meta", type=Path, required=True,
             help="蛋白元数据文件路径 (TSV格式)")
@click.option("--epitope-ratios", type=Path, required=True,
             help="表位比率文件路径 (CSV格式)")
@click.option("--output", "-o", type=Path, default="enhanced_merged_data.tsv",
             help="整合数据输出路径 (TSV格式)")
def merge_data(pdb_mapping: Path, proteins_meta: Path, epitope_ratios: Path, output: Path):
    """整合PDB映射、蛋白元数据和表位比率数据（保留所有表位数据）"""
    merger = EnhancedDataMerger(
        pdb_mapping=pdb_mapping,
        proteins_meta=proteins_meta,
        epitope_ratios=epitope_ratios,
        output=output
    )
    
    try:
        merger.validate()
        click.secho("🔄 开始整合数据...", fg='green')
        merger.merge_data()
        click.secho(f"✅ 整合完成！结果保存至 {output}", fg='green', bold=True)
    except Exception as e:
        click.secho(f"❌ 错误: {str(e)}", fg='red')
        
# ================= 增强版TreeMap可视化命令 =================
@cli.command()
@click.option("--input", "-i", type=Path, required=True,
             help="整合后的TSV数据文件路径")
@click.option("--output", "-o", type=Path, default="antigen_treemap.png",
             help="输出文件路径（默认为PNG格式）")
@click.option("--species", type=str, default=None,
             help="按物种过滤（如'Homo sapiens'）")
@click.option("--min-ratio", type=float, default=0.0,
             help="最小表位比率阈值（0-1之间）")
@click.option("--format", "output_format", type=click.Choice(['png', 'svg', 'jpg', 'pdf']),
             default='png', help="输出文件格式（默认为PNG）")
def plot_treemap(input: Path, output: Path, species: str, min_ratio: float, output_format: str):
    """生成抗原表位比率的TreeMap可视化（增强版）"""
    visualizer = TreeMapVisualizer(
        input_file=input,
        output_file=output,
        species_filter=species,
        min_ratio=min_ratio,
        output_format=output_format
    )
    
    try:
        click.secho("🌳 正在生成TreeMap可视化...", fg='green')
        visualizer.execute()
        click.secho(f"✅ 可视化已成功保存至: {output}", fg='green', bold=True)
    except Exception as e:
        click.secho(f"❌ 生成TreeMap失败: {str(e)}", fg='red')

if __name__ == "__main__":
    cli()
