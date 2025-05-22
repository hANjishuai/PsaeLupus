"""
抗原分析工具集 - 整合版

子命令列表
calculate-ratios : 分块计算表位比率

"""

import csv
import time
import click
import functools
from pathlib import Path
import matplotlib.pyplot as plt
from typing import Dict, Generator, Tuple
from collections import defaultdict

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


if __name__ == "__main__":
    cli()
