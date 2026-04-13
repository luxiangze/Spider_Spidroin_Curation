# 001 号物种集成测试 — *Allagelena difficilis*

## 概述

本测试套件以人工鉴定结果为 ground truth，验证蜘蛛丝蛋白自动鉴定流程（typing_agent）
对 001 号物种（*Allagelena difficilis*，缩写 `Aldi`）的运行正确性与结果准确性。

所有测试直接调用 `_process_species()` 函数（绕过 CLI），并关闭基因组注释重叠、
密码子检查和 RNA-seq 证据（`skip_anno`、`skip_codon`、`skip_rna`），
以避免外部路径依赖并将运行时间控制在约 30 秒内。

---

## Ground Truth

`ground_truth.csv` 包含 23 条 spidroin 基因位点记录，由人工在 JBrowse2 中
逐条审查（2026-03-22）。每条记录附有**置信度评分**（1–5）：

| 评分 | 含义 |
|------|------|
| 5 | 无歧义 — NTD 和 CTD 命中质量高，有强烈 miniprot 支持 |
| 4 | 可信 — NTD+CTD 配对明确 |
| 3 | 可能正确 — 证据较弱 |
| 2 | 部分 — 仅检测到 CTD（约 330 bp） |
| 1 | 极低置信度 |

截图证据存储在 `attachments/` 文件夹中。

**结果汇总：**

| 类型 | 完整性 | 数量 | 评分范围 |
|------|--------|------|---------|
| AcSp | Full_length | 5 | 3–4 |
| AcSp | C-terminal  | 7 | 1–2 |
| PySp | Full_length | 3 | 5 |
| CySp | Full_length | 2 | 4 |
| MiSp | Full_length | 5 | 5 |

---

## 关键测试案例：chr07 复杂窗口

区域 `chr07:71,243,836–71,465,633` 经抑制处理后包含 **1 个 NTD + 6 个 CTD** 命中，
即 `agents/samples/sample7.png` 中展示的"多个相邻 CTD"典型场景。

人工审查结论：
- **spid00004–00008**：5 个 C-terminal（各约 330 bp，评分 2，非全长）
- **spid00009**：AcSp 全长位点 `chr07:71,442,768–71,465,633`（评分 4）  
  — NTD 位于 71,465,162，与**最近的 CTD**（71,442,774，跨度 22,860 bp）正确配对，
  该区间内有 11 个非重叠 miniprot 比对作为强证据支持。

本测试案例专门验证配对算法是否能避免将 NTD 与更远的 CTD 错误配对。

---

## 运行方式

```bash
# 快速测试（不调用 LLM，约 30 秒）
pixi run pytest tests/test_species_001/ -v

# 含 LLM 配对测试（需 .env 中有 ANTHROPIC_API_KEY，约 2 分钟）
pixi run pytest tests/test_species_001/ -v --run-llm
```

---

## 输出文件

测试运行后，`output/` 文件夹包含：

| 文件 | 说明 |
|------|------|
| `001.no_llm.tsv` | 使用贪心算法的结果（无 LLM） |
| `001.llm.tsv` | 使用 LangGraph + Claude 配对的结果（仅 --run-llm 时生成） |

输出列包含标准飞书字段及 `llm_confidence`、`llm_reasoning`，可直接打开检查 LLM 决策。

---

## 测试用例说明

| 测试函数 | 说明 |
|---------|------|
| `test_pipeline_runs` | 冒烟测试 — 流程无错误完成，所有位点结构合法 |
| `test_high_confidence_full_length` | 人工标注评分 ≥ 4 的 13 个全长位点必须全部被检出 |
| `test_ambiguous_window_chr07` | chr07 复杂窗口正确配对；其余 ≥ 4 个 CTD 保留为 C-terminal |
| `test_no_false_full_length_in_ambiguous_window` | NTD 不得与更远的错误 CTD 配对为全长 |
| `test_llm_pairing` *（需 --run-llm）* | LLM 路径给出 `high` 置信度，且包含 spanning miniprot 证据标签 |
| `test_output_tsv_written` | 输出 TSV 文件存在且数据行数 ≥ 15 |
