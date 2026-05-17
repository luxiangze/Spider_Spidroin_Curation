# 蜘蛛丝蛋白自动鉴定流程：方法介绍与初步结果

**汇报人：** 郭永康
**日期：** 2026年4月8日  
**项目：** 蜘蛛 Silkome 蜘蛛丝蛋白自动注释流程

---

## 一、研究背景与目标

蜘蛛丝蛋白（spidroin）是蜘蛛丝腺分泌的主要功能蛋白，不同类型（MaSp、MiSp、Flag 等）对应不同的丝腺和力学性能。在大规模蜘蛛基因组注释中，spidroin 基因的准确鉴定与分型是构建 silkome 数据库的核心步骤。

**现有工作流的局限性：** 目前的流程依赖人工在 JBrowse2 中逐一审查 nHMMER 结构域搜索结果和 miniprot 参考蛋白比对结果，再手动录入飞书多维表格。该流程在处理 100 余个物种的基因组时效率极低，且人工判断标准难以统一。

**本工作的目标：** 构建一套自动化的 spidroin 鉴定与分型流程（`typing_agent`），实现对 nHMMER 和 miniprot 结果的整合分析，输出结构化的 TSV 结果文件，同时对低置信度结果打标记以供后续人工复核，从而替代绝大部分手动审查工作。

---

## 二、流程设计与方法

### 2.1 整体架构

流程主要分为五个阶段，如下图所示：

```
nHMMER GFF（NTD/CTD 结构域）
        ↓
   位点组装（Locus Assembly）
        ↓
miniprot GFF（参考蛋白比对）  ─→  证据整合（Evidence Integration）
基因组注释（genome annotation） ─↗
RNA-seq 覆盖度（BigWig）      ─↗
        ↓
   分型与完整性评估（Typing & Completeness）
        ↓
   置信度计算 → 结果输出（TSV / GFF / FASTA）
```

### 2.2 主要方法步骤

**（1）nHMMER 结果解析**

nHMMER 针对核苷酸序列进行 HMM profile 搜索，每种 spidroin 类型均构建了 NTD 和 CTD 两套 profile（命名规则：`{类型}_{结构域}.hmm`，如 `MaSp1_NTD.hmm`）。解析后得到每个结构域命中的位置、链方向和 E-value。

筛选阈值：
- E-value < 1×10⁻⁴（初步过滤）
- E-value < 1×10⁻¹⁰（高置信度阈值）
- score > bias（排除低复杂度假阳性）

**（2）位点组装（Locus Assembly）**

采用贪心配对策略，在同一条 scaffold 上将 NTD 和 CTD 命中进行两两匹配，组装候选 spidroin 基因位点。配对约束条件：

- NTD 与 CTD 中点距离：1 kb ≤ span ≤ 150 kb
- 同侧链方向
- 优先保留 E-value 最优命中，抑制重叠命中
- 泛用型 profile（`Spidroin_NTD/CTD`）与特异型命中重叠时被抑制

未能成功配对的命中单独保留，标注为 N-terminal-only 或 C-terminal-only 候选。

**（3）多证据整合**

- **miniprot 比对：** 找到与候选位点重叠的最优参考蛋白比对（优先 rank=1 且 positive score 最高），记录参考蛋白名称、identity 和 positive 分值。
- **基因组注释：** 查询已有注释基因（gene_id）与候选位点的重叠比例。
- **RNA-seq：** 从 BGI 短读长和 Oxford Nanopore 长读长 RNA-seq BigWig 文件中提取候选区域的转录信号强度。
- **密码子验证：** 检查 NTD 起始密码子和 CTD 终止密码子是否存在。

**（4）类型分型（Type Classification）**

分型规则从外部 YAML 文件（`docs/typing_rules.yaml`）加载，不硬编码于代码中，方便后续调整。分型逻辑如下：

- NTD profile 和 CTD profile 均映射到标准类型名称（共 16 种主要类型）
- 若两侧类型一致 → 直接确定类型
- 若存在特异性层级关系（如 `MaSp1` > `MaSp`）→ 取更特异的类型
- 若属于同一家族的不同亚型（如 `MaSp1` vs `MaSp2`）→ 退化为家族名（`MaSp`），标记 `needs_review=True`
- 若属于不同家族 → 记录冲突（如 `MaSp/MiSp`），标记 `needs_review=True`

目前支持的 spidroin 类型：MaSp、MaSp1、MaSp2、MaSp2B、MaSp3、MaSp3B、MiSp、PySp、AgSp、AgSp1、AgSp2、Flag、Pflag、TuSp、AcSp、CySp、CrSp

**（5）完整性评估（Completeness Assessment）**

| 类别 | 判定条件 |
|------|---------|
| Full_length | NTD 和 CTD 均存在，且均达到高置信度 E-value 阈值 |
| N-terminal | 仅 NTD 存在（高置信度） |
| C-terminal | 仅 CTD 存在（高置信度） |
| Repeat_only | 无高置信度结构域，但有 miniprot 比对支持 |

**（6）置信度评级**

综合所有证据进行评级（高 / 中 / 低）：
- **High：** NTD 和 CTD 均高置信度，类型无冲突
- **Medium：** 类型存在冲突，或仅一侧达到高置信度
- **Low：** 仅有弱证据或仅有注释/RNA-seq 支持

### 2.3 输出格式

每个物种生成四类文件：

| 文件 | 内容 |
|------|------|
| `{species}.tsv` | 主结果表（32 列，含飞书字段和诊断字段） |
| `{species}.gff` | GFF3 格式的 spidroin 基因特征文件 |
| `spidroin_full_length.fasta` | 全长候选的基因组序列（FASTA） |
| `hints.gff` | Augustus 基因预测提示文件（起止密码子标记） |

---

## 三、当前处理结果

### 3.1 总体概况

目前已完成 **127 个蜘蛛物种** 的自动化分型，共鉴定候选 spidroin 位点 **2157 个**。

### 3.2 完整性分布

| 类别 | 数量 | 占比 |
|------|------|------|
| Full_length | 1743 | 80.8% |
| N-terminal only | 235 | 10.9% |
| C-terminal only | 179 | 8.3% |
| Repeat_only | 0 | 0.0% |
| **合计** | **2157** | **100%** |

全长候选占比超过 80%，说明多数 spidroin 基因可被完整鉴定。目前尚未检出 Repeat_only 类型，主要原因是该类别的判定逻辑依赖 repeat 区域检测，尚未实现（见第四节）。

### 3.3 置信度分布

| 置信度 | 数量 | 占比 |
|-------|------|------|
| High | 1292 | 59.9% |
| Medium | 865 | 40.1% |
| **合计** | **2157** | **100%** |

约 60% 的候选位点可自动确认（高置信度，无需人工复核），相较于原有的全人工审查流程，可大幅减少工作量。其余 40.1%（865 个）标记为 `needs_review=True`，主要由类型冲突或单侧置信度不足引起。

### 3.4 类型分布（Top 15）

| Spidroin 类型 | 数量 | 备注 |
|--------------|------|------|
| MiSp | 461 | 管状腺丝蛋白，最丰富 |
| AcSp | 285 | 管壶腺丝蛋白 |
| MaSp1 | 262 | 壶状腺丝蛋白 1 型 |
| MaSp | 262 | 壶状腺泛型（NTD/CTD 特异性不同） |
| Flag | 134 | 鞭状腺丝蛋白 |
| PySp | 125 | 梨状腺丝蛋白 |
| CySp | 107 | 葡萄状腺丝蛋白 |
| MaSp/MiSp | 102 | 类型冲突（NTD vs CTD 归属不同家族） |
| MaSp2 | 81 | 壶状腺丝蛋白 2 型 |
| AgSp2 | 80 | 集合腺丝蛋白 2 型 |
| MaSp2B | 60 | — |
| AgSp1 | 53 | — |
| MaSp/AcSp | 33 | 类型冲突 |
| MaSp3B | 23 | — |
| MaSp3 | 23 | — |

类型冲突条目（如 `MaSp/MiSp`，`MaSp/AcSp` 等）合计约 156 个（7.2%），均已标记 `needs_review=True`，需结合人工判断进行最终归类。

### 3.5 密码子边界验证

TSV 中包含 `ntd_start_codon` 和 `ctd_stop_codon` 两列，记录 NTD 起始密码子与 CTD 终止密码子是否匹配。全长候选（Full_length）中密码子边界完整的条目可进一步提升置信度。

---

## 四、存在的问题与待改进之处

### 4.1 Repeat_only 类型未检出

当前流程缺少 repeat 区域检测模块。对于仅含重复单元、不含 NTD/CTD 的 spidroin 截段，目前无法识别。后续需引入基于 miniprot 比对结果中 repeat 区域的判断逻辑。

### 4.2 miniprot 置信度阈值待确认

当前 miniprot 的高置信度阈值（positive score ≥ 0.6；identity ≥ 70%）为初步设定值，需结合人工审查结果进行验证与调整。

### 4.3 跨家族类型冲突（156 个条目）

部分候选位点的 NTD 和 CTD 分别命中不同科的 spidroin profile（如 MaSp 和 MiSp），可能反映真实的基因融合或 domain shuffling 现象，也可能是 HMM profile 交叉匹配的假阳性。需人工逐一核查。

### 4.4 泛用型 profile 的 "污染"

泛用 `Spidroin_NTD/CTD` profile 在部分物种中出现大量命中，导致部分条目无法精确分型（归为 `unknown` 或 `needs_review`）。后续可考虑针对这些物种的 spidroin 序列重新构建特异性 profile。

---

## 五、后续工作计划

1. **人工抽样复核（近期）**：随机抽取约 50 个 `needs_review=True` 的条目，在 JBrowse2 中人工审查，评估当前分型准确率，并据此调整阈值参数。

2. **Repeat_only 检测（近期）**：实现基于 miniprot target 区域的重复单元覆盖度分析，补全 Repeat_only 类型的识别逻辑。

3. **飞书 API 对接（中期）**：在现有 TSV 输出基础上，开发飞书多维表格写入模块，实现高置信度结果的自动录入。

4. **跨物种统计分析（中期）**：在全部 134 个物种完成鉴定后，对各 spidroin 类型的分布进行跨物种比较分析，探讨 silkome 多样性的演化规律。

5. **流程验证（持续）**：以已知物种（文献有详细注释）的结果为标准，对流程的灵敏度和精确率进行定量评估。

---

## 附录：关键参数汇总

| 参数 | 值 |
|------|----|
| nHMMER E-value 初步筛选 | < 1×10⁻⁴ |
| nHMMER 高置信度阈值 | < 1×10⁻¹⁰ |
| miniprot 高置信度 positive score | ≥ 0.6 |
| 位点 span 范围 | 1 kb – 150 kb |
| 支持的 spidroin 类型数 | 17 种 |
| 已处理物种数 | 127 |
| 已鉴定候选位点数 | 2157 |
| 高置信度自动确认率 | 59.9% |
| 需人工复核率 | 40.1% |
