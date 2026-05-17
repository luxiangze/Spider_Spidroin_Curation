# 蛛丝蛋白数据表 - 列名清单

| # | 列名 | 类型 | 说明 |
|---|------|------|------|
| 1 | Spidroin_ID | 文本 | 主键，蛛丝蛋白 ID |
| 2 | Chr | 文本 | 染色体 |
| 3 | Start | 数字 | 起始位置 |
| 4 | End | 数字 | 终止位置 |
| 5 | Stran | 单选 | 链方向（+ / -） |
| 6 | Spidroin_type | 单选 | 蛛丝蛋白类型 |
| 7 | Full_length | 单选 | 是否全长（True / False） |
| 8 | Length | 公式 | 自动计算长度（End - Start + 1） |
| 9 | Hint_type | 单选 | 提示类型（Full_length, C-terminal, N-terminal, Internal） |
| 10 | Species | 单选 | 物种来源 |
| 11 | Note | 文本 | 备注 |
| 12 | attachment | 附件 | 附件字段 |
| 13 | Scoring | 评分 | 评分（1-5 分） |
| 14 | People | 修改人 | 自动记录最后修改人 |
| 15 | Update_time | 修改时间 | 自动记录最后修改时间 |

## 字段详情

### 单选字段选项

**Spidroin_ID 格式：**`{species}_spid_{number}`，例如：`Arve_spid_00001`

**Stran（链方向）:**
- `+`
- `-`

**Spidroin_type（蛛丝蛋白类型）:**
Flag, AgSp2, PySp, MaSp1, MaSp, MiSp, MaSp2, MaSp3, AcSp, CySp, MaSp3B, AgSp1, Unkown, MaSp2B, CrSp, hypo, flag, Pflag

**Full_length（是否全长）:**
- True
- False

**Hint_type（提示类型）:**
- Full_length
- C-terminal
- N-terminal
- Internal

**Species（部分物种示例）:**
013.Evarcha sp, 017.Hippasa lycosina, 031.Pandercetes sp, 034.Pholcus sp, 045.Scorpiops zhui, 049.Songthela sp, 064.Araneus ventricosus, 079.Heteropoda venatoria, 106.Pardosa pseudoannulata, 119.Trichonephila clavata 等...