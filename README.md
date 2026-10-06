# 生物序列分析平台

这是一个只使用 Python 标准库的离线平台。第一版保存参考序列、输入数据集（含非法记录）、参数和运行批次，并为后续变异分析输出坐标、CIGAR、得分和变异事件。

## 能力

- 本地 Smith-Waterman（`--mode local`，默认）和 Needleman-Wunsch（`--mode global`）。
- 支持替换、删除、插入和仿射空位惩罚。
- 回溯所有等分最优候选；重复区域不会只报告第一个位置。
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

## 环状拓扑

- 环状分子模式（`--topology circular`）：枚举参考序列的每个旋转起点并分别做线性比对，跨越记录末尾的匹配是一段完整比对而非碎片；方向（正/反向互补）与环上起点（`rotation_offset`）都会报告，坐标始终回映到原始输入坐标。
- 环状模式只把查询序列使用一次：得分来自“旋转后的参考 × 原始查询”的单次线性比对，短序列或端部缺口不会被重复计入以虚增得分。
- 重复片段产生的多个同分起点会全部保留为并列候选（按 `strand`、`rotation_offset`、终点坐标、路径排序）。
- `wraps_reference` 标记坐标是否跨过记录末尾；跨界的替换/插入缺失在 `variants` 中以原始输入坐标给出。
- 同一份输入可分别以 `--topology linear` 与 `--topology circular` 运行，得到不同 `run_id`，便于比较两种解释。
- 人工选定候选：`bioseq select <project> <run_id> <record_id> <index>` 将 `selected_alignment_index` 写回该批次；重新读取后候选顺序、方向和选定起点保持不变。

## 功能区注释与投影

- 用 `bioseq feature-add <project> <reference_id> --name N --start S --end E [--strand +/-]` 在原始参考上标记有方向区间；`start > end` 表示跨过记录起点，该区间始终作为一个整体处理，不会拆成两段。
- `bioseq project <project> <run_id> <record_id> [--index I]` 在每个并列候选解释下把功能区投影到查询序列：输出 `query_start`/`query_end`（原始查询坐标，`query_wraps_origin` 标记跨界）、`query_strand`、覆盖的参考位置和落在功能区内的变异。所有同分解释各自的投影全部保留，不自动挑选影响最小的一种。
- 变异归属规则：替换按参考位置判定；缺失与功能区有交集即报告，`covers_feature` 标记是否覆盖整个功能区——覆盖时功能区状态为 `deleted`，不再有对应查询范围；插入 `after_position=p` 仅当 p 与 p+1（环状意义下）都在区间内时记为 `internal`，只沾一边记为 `boundary`，因此边界插入不会在相邻两个功能区里都被算成完整内部变异。
- 反向互补解释下，投影坐标、方向和区间先后关系都回映到原始输入序列（`query_strand=reverse`，坐标为原始查询坐标）。
- 改选解释前可用 `bioseq select ... --preview` 查看哪些功能区的投影需要重算（`features_to_recompute`）；正式 `select` 的输出也带 `features_recomputed`。功能区注释按 `reference_id` 保存在 `features/` 下，与批次结果相互独立，改选后旧投影仍可按候选索引重现比较。

## 参考版本链与特征追踪

- `bioseq chain-create <project> <reference_id> [--chain-id ID]` 以原始参考为第 0 版建立连续版本链，保存在 `chains/` 下。
- `bioseq chain-edit <project> <chain_id> rotate|delete|insert ...` 追加一个新版本：先旋转（`--offset`）再删除（`--start/--end`）或插入（`--after/--bases`），每个版本记录逐位置来源（carried 指向上一版本坐标，inserted 标记产生它的编辑版本与偏移）。
- `bioseq chain-trace <project> <chain_id> <feature_id> [--version K]` 把原始参考上的有方向特征（允许 `start > end` 跨起点）追到新链的每个中间版本：输出保留原特征定义与分段顺序（`segments` 按特征方向排列，含每段来源），不是排序位置集合；部分删除标记 `partially_deleted` 并列出 `deleted_sources`/`deleted_source_positions`，整体删除明确标记 `deleted`；插入编辑附带 `insertion_impact`（internal/boundary/outside，沿用边界插入归属规则）。`--version K` 只导出到第 K 版的中间映射。
- 已归属特征的内部插入片段会成为特征成员并在后续版本继续参与映射：每个成员的来源是原子描述符 origin（原始位置来源，字段 position）或 inserted（插入版本来源，字段 edit_version/offset），见每个版本快照的 sources 与各分段的 sources/source_positions（插入碱基在后者记 null）。后续旋转保持成员的方向顺序，后续删除若落在插入片段上会在 deleted_sources 里列出对应的插入版本来源。
- 边界插入（只沾特征一侧）与外部插入不归入特征：insertion_impact.attributed=false、inserted_sources 为空，之后的旋转或删除也不会把这些碱基带进特征。
- `bioseq chain-undo <project> <chain_id>` 撤回末次编辑（丢弃最后一个版本）并写回。
- `bioseq chain-verify <project> <chain_id>` 用直接映射（编辑变换一次性组合）与逐版本来源组合核对每个版本的全部原子（原始位置 + 各版本插入的碱基），不一致时列出差异。
- `bioseq chain-plan <project> <chain_id> --candidates edits.json [--features F1,F2] [--allowed-breaks breaks.json] [--confirm-chain-id ID]` 在环形原始坐标上研究多个旋转、删除和插入候选的可共同应用组合。候选必须有稳定 `id`；两个候选即使插入完全相同的碱基，`id`/offset 也代表不同插入来源。删除用 `start/end`（允许跨界）、`positions` 或 `deleted_sources` 指向 origin/inserted 原子；插入可用 `after_source` 指向插入版本来源。
- 规划器枚举满足受保护有向特征及允许破坏范围的组合，优先保留最多候选，再比较删除和插入碱基数最少者；结果包含全部可行组合、共享删除范围的候选、连续插入位点，以及冲突的实际来源与候选。冲突判定只基于原子来源、旋转后的原始位置、特征方向和跨起点分段顺序，不依赖候选登记顺序。
- 传入 `--confirm-chain-id` 后，最佳组合会从第 0 版重新生成独立版本链并保存；该链可继续用 `chain-trace` 检查反向/跨界特征归属，用 `chain-verify` 检查共享删除、插入来源和连续插入后的最终序列。
- 保存重开后段来源与插入归属保持不变，trace 与 verify 结果与重开前一致。
