# 生物序列分析平台

这是一个只使用 Python 标准库的离线平台。第一版保存参考序列、输入数据集（含非法记录）、参数和运行批次，并为后续变异分析输出坐标、CIGAR、得分和变异事件。

## 能力

- 本地 Smith-Waterman（`--mode local`，默认）和 Needleman-Wunsch（`--mode global`）。
- 支持替换、删除、插入和仿射空位惩罚。
- 回溯所有等分最优候选；重复区域不会只报告第一个位置。
- 环状分子模式（`--topology circular`）：对环状查询枚举所有旋转起点（而非加倍序列），
  同时搜索正向与反向互补，结果给出 `strand` 与 `query_rotation`（环上起点），
  跨过记录末尾的匹配是一条完整比对而不是割裂片段；坐标始终以原始输入给出（可环绕）。
  每个碱基在一次比对中至多使用一次，短序列或端部缺口不会被重复计入得分。
- 环状模式下同分旋转起点全部保留，可用 `bioseq select <项目> <run_id> <记录id> <候选序号>`
  人工选定候选，选定结果写回批次并在重新读取后保持。
- 环状候选的变异事件附带原始查询坐标（`query_position`/`query_start`/`query_end`/
  `query_after_position`），跨界替换与插入缺失可直接定位到原始记录。
- 非法输入按记录隔离，原因写入同一批次结果，不影响其他序列。
- 数据集、参考序列、参数和运行结果分别保存为内容寻址 JSON。
- 同一项目、输入和参数重新运行会得到同一个确定性 `run_id` 并复用已保存结果。
- 不调用在线生物数据库。

## 快速使用

```bash
python -m venv .venv
. .venv/bin/activate
pip install -e .

bioseq init demo-project
bioseq align demo-project \
  --reference examples/reference.fasta \
  --sequences examples/sequences.fasta \
  --type dna --mode global

bioseq show-run demo-project run-<返回的批次哈希>
```

也可以不安装包，直接设置 `PYTHONPATH=src` 后运行：

```bash
PYTHONPATH=src python -m bioseq_platform.cli init demo-project
```

## 存储结构

- `references/`：规范化参考序列和元数据。
- `datasets/`：有效与非法输入记录及非法原因。
- `parameters/`：比对模式、得分、空位惩罚和候选上限。
- `runs/<run_id>/manifest.json`：该批次的完整追溯信息。
- `runs/<run_id>/results.json`：每条记录的状态、候选比对、CIGAR 和变异。

坐标为 1 基闭区间；CIGAR 使用 `=` 匹配、`X` 替换、`D` 参考删除、`I` 查询插入。

每个候选的身份是其在原始序列上的坐标路径（`reference_coordinates`/`query_coordinates`，gap 列记 0），而非对齐后的文字，因此重复区域中文字相同、来源不同的并列解都会保留。局部与全局比对都会同时搜索正向与反向互补链，候选带 `strand`（`forward`/`reverse`）标记；反向候选的 `query_start`/`query_end` 与查询坐标按原始（未反向）序列给出。正反两个方向得分相同的并列候选会全部保留（按 `strand`、终点坐标、路径排序），不会静默丢弃任一方向；候选的 `score` 即为方向取舍的依据。蛋白质只搜索正向。
