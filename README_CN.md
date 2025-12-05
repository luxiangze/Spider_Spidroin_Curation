# Spider Spidroin Curation

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

**[English](README.md)**

蜘蛛丝蛋白（Spidroin）自动化注释与表达分析流程。

## 快速开始

```bash
# 安装依赖
pixi install

# 激活环境
pixi shell

# 运行 Spidroin 注释流程（Jupyter Notebook，需要先设置配置参数）
jupyter lab notebooks/Automated_spidroin_annotation.ipynb

# 运行 BGI 短读长 RNA-seq 比对流程
cd workflow/BGI_RNA_align
pixi run snakemake --cores 80 --resources mem_gb=240

# 运行 ONT 长读长 RNA-seq 转录本组装（需要先完成 BGI 流程）
cd workflow/ONT_RNA_align
pixi run snakemake --cores 80 --resources mem_gb=240
```

## 主要功能

### Spidroin 自动化注释
- ✅ nhmmer 搜索 Spidroin N/C 端序列
- ✅ 多物种批量分析（10 个蜘蛛基因组）
- ✅ Augustus 基因结构预测
- ✅ 蛋白序列提取

### RNA-seq 分析
- ✅ BGI 短读长 STAR 比对流程 (`workflow/BGI_RNA_align`)
- ✅ ONT 长读长 FLAIR 转录本组装 (`workflow/ONT_RNA_align`)
- ✅ FastQC + fastp 质控
- ✅ BigWig 可视化文件生成
- ✅ MultiQC 汇总报告

## 项目结构

```
spider_silkome/
├── data/
│   ├── raw/                       # 原始数据
│   │   ├── spider_genome/         # 蜘蛛基因组 (fa + gff)
│   │   ├── BGI_RNA_10samples/     # BGI 短读长 RNA-seq 数据
│   │   └── ONT_DRS_10samples/     # ONT 长读长 RNA-seq 数据
│   ├── interim/                   # 中间处理结果
│   └── processed/                 # 最终输出
│
├── workflow/
│   ├── BGI_RNA_align/             # BGI 短读长 STAR 比对流程
│   │   ├── Snakefile
│   │   ├── config/
│   │   └── rules/
│   └── ONT_RNA_align/             # ONT 长读长 FLAIR 转录本组装
│       ├── Snakefile
│       ├── config/
│       └── rules/
│
├── spider_silkome_module/         # Python 模块
│   ├── config.py                  # 路径配置
│   ├── features.py                # 工具函数
│   └── plots.py                   # 可视化
│
├── notebooks/
│   └── Automated_spidroin_annotation.ipynb  # Spidroin 注释主流程
├── scripts/
│   └── analyse_spidroins.py       # Spidroin 分析脚本
└── Makefile
```

## 依赖管理

项目使用 [pixi](https://pixi.sh/) 管理依赖，主要工具包括：

- **Spidroin 注释**: nhmmer (HMMER), Augustus, BioPython
- **短读长比对**: STAR, samtools
- **长读长比对**: minimap2, FLAIR
- **质控**: FastQC, fastp, MultiQC
- **可视化**: deeptools (bamCoverage)
- **流程**: Snakemake, Jupyter

## 许可证

- **许可证**: 详见 [LICENSE](LICENSE) 文件
- **维护者**: 郭永康
- **模板**: 基于 [Cookiecutter Data Science](https://cookiecutter-data-science.drivendata.org/)
