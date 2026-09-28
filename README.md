# 生物序列分析平台

纯 Python 标准库实现的离线小规模序列分析平台。首个能力是 Smith-Waterman 局部序列比对，使用线性缺口罚分，完整保留最佳得分下的并列候选。

## 数据留存

默认数据目录为 `.bioseq/`：

- `references/`：规范化参考序列和内容哈希。
- `datasets/`：输入数据集和内容哈希。
- `parameters/`：匹配、错配和缺口参数。
- `batches/<batch_id>/`：批次清单、有效结果和隔离记录。

批次 ID 由算法版本、数据集 ID、参考 ID、参数 ID 和候选上限的稳定哈希生成；不写入时间戳。重新运行同一项目会生成同一个批次目录和相同结果，因此可用于后续变异分析追溯。

## 使用

```bash
PYTHONPATH=src python -m bioseq.cli \
  --store .bioseq \
  align \
  --reference examples/reference.json \
  --dataset examples/dataset.json
```

查看批次：

```bash
PYTHONPATH=src python -m bioseq.cli --store .bioseq show <batch_id>
```

支持 DNA、RNA 和蛋白质字符集。数据集可使用 JSON，或包含 `id,sequence,sequence_type` 列的 CSV。无效或重复 ID 的记录写入 `quarantine.json`，不会中断同批其他记录。

## 结果坐标

坐标是零基、左闭右开区间。`aligned_query` 中的 `-` 表示查询相对参考发生缺失（CIGAR `D`）；`aligned_reference` 中的 `-` 表示查询插入（CIGAR `I`）。重复区域中的每个同分位点都会作为候选保留，不会只返回第一个位置。

## 验证

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```
