# 推理文件与区域选择的字段契约

一个样本包含一张实验图、一种受体条件、一个去重的配体身份列表。下文用 **I** 表示规范 SMILES 身份数，**B** 表示 Find 输出的可考察区域数。一个身份可以在多个位置出现。样本与区域联合用 `(sample_id, blob_id)` 标识；样本内的真实拷贝用 `occurrence_id` 区分。命令和可运行例子见[使用说明](README.md)。

以下是同一个输出目录的契约，文件随完成的阶段写出。完成阶段的区域和姿态均可从这些文件定位。

```text
## <科学产物>
reviewed/
  find.json                    # 各区域的原编号、中心和定位分数
  match.json                   # 同一区域的前景概率与身份得分
  match_model.json             # 模型的前景阈值与可编码身份下标，用于重建准确率分母
  find/inputs/                 # 同网格的实验图、模拟图和受体数组，评估在该网格上重建参考
  find/probability/            # 整图概率和网格几何，供语义 PRAUC 评价
  find/blobs/F2_blobs.npz       # F2 语义阈值提取的完整区域，供独立空间归属
  find/centered/F2_centered.npz # 局部特征、评分和选择状态，供 Match 和定位评价
  plan.json                    # 每个区域的计划与依据
  build.json                   # 实际提交的请求和姿态位置
  build/17/poses.sdf           # 区域 17 的成功姿态，17 代表一个实际提交的 blob_id
  build/17/candidates.json     # 各原采样编号与成功 SDF 位置
  build/17/ranking.json        # 逐姿态综合评分与首选编号
  build/17/confidence.npz      # 与成功姿态顺序对应的置信度
## <其他文件>
reviewed/
  run.json                     # 输入路径与本次命令参数（运行诊断）
```

## 样本、身份与实例

### `sample.json`

| 字段 | 类型或形状 | 科学含义 | 构造示例 |
| --- | --- | --- | --- |
| `sample_id` | 字符串 | 唯一标识样本及受体条件，阶段间保持相同 | `"demo-ca2"` |
| `experimental_map` | 路径字符串 | 实验密度图，体素数组按 ZYX 编号 | `"data/experiment.mrc"` |
| `receptor_source` | 字符串 | 选择 Find 的冻结校准。已有沉积受体填 `deposited`，有序列或无序列的 CryoAtom2 受体填 `ca2`。旧输入缺省为 `ca2` | `"ca2"` |
| `receptor_cif` | 路径字符串 | 与实验图共用世界 XYZ 坐标系的受体 mmCIF。也可提供完整复合物，Find 从中提取受体原子。 | `"data/ca2.cif"` |
| `simulated_map` | 路径字符串 | 用同一受体模拟的密度，与实验图使用相同的原点、间距和网格 | `"data/ca2_density.mrc"` |
| `map_convention` | 可选字符串 | `mrc` 用同一套实空间坐标约定读取两张 MRC；`adaligand` 复用原训练数据的读取方式。缺省时，模拟图为 AdaLigand NPY/NPZ 资产取后者，成对 MRC 取前者。细节见 [Find 网格说明](ligandseek/find/README.md) | `"mrc"` |
| `identities` | 对象列表，长度 I | 首次出现顺序的不同规范 SMILES，每个身份占一格 | `[{"smiles":"CCO","type_tag":"small_molecule"}]` |
| `identities[].smiles` | 字符串 | 在本样本内唯一的身份 | `"CCO"` |
| `identities[].ccd_id` | 可选字符串 | 用户给定的通用 CCD 组分名，可辅助识别糖、核酸和肽类 | `"NAG"` |
| `identities[].type_tag` | 可选字符串 | 已知类别；未知时缺省或为 null | `"small_molecule"`、`"sugar"` |
| `identities[].embedding` | 可选路径字符串 | SMI-TED 的 float32 (768,) 向量，NPY 直接存向量，NPZ 从 `embedding` 读取。仅制定计划时可省略 | `"data/ethanol.npy"` |

统一身份调用 RDKit 的 canonical non-isomeric SMILES。四面体标记、双键方向、同位素标号均不进入身份键。构造例 `C[C@H](O)C(=O)O` 和 `C[C@@H](O)C(=O)O` 共用一个 `CC(O)C(=O)O`，Match 产生一个分数，自动和手工 Build 都以该键构图。姿态仍保存模型生成的三维坐标。不能解析的原文映射为 `null`。

从带有沉积配体的 mmCIF 准备输入时，先提取仅含受体的原子并用它们生成模拟图。完整沉积结构可用于准备独立的参考实例。Find、Match 和 Build 使用同一受体。样本相对路径从 `sample.json` 所在目录解析。

### `reference.json`

窄范围可读取已判定的实例归属。一个完整构造例子如下。

```json
{"sample_id":"demo-ca2","identity_types":{"CCO":"small_molecule"},"occurrences":[{"occurrence_id":0,"smiles":"CCO","type_tag":"small_molecule","eligible":null,"complete_heavy_atoms":false}],"assignments":[{"blob_id":17,"occurrence_id":0}]}
```

| 字段 | 类型或形状 | 含义 | 构造示例 |
| --- | --- | --- | --- |
| `sample_id` | 字符串 | 样本及受体条件标识 | `"demo-ca2"` |
| `identity_types` | 可选字典，键为 SMILES，值为类别 | 按身份补充 `sample.json` 中缺少的类别 | `{"CCO":"small_molecule"}` |
| `occurrences` | 长度 G 的列表 | G 为样本中的真实配体拷贝数，每项描述一个拷贝 | 见上例，G = 1 |
| `occurrences[].occurrence_id` | 非负整数 | 该真实拷贝在样本内的唯一编号 | `0` |
| `occurrences[].smiles` | 字符串或 null | 该拷贝的精确 SMILES | `"CCO"` |
| `occurrences[].type_tag` | 类别字符串或 null | 该拷贝的配体类别 | `"small_molecule"` |
| `occurrences[].eligible` | 可选布尔值或 null | 该拷贝是否满足论文 Stage3 总资格 | `null` |
| `occurrences[].single_residue_ccd` | 可选布尔值或 null | 是否为单残基 CCD 组分 | `true` |
| `occurrences[].noncovalent` | 可选布尔值或 null | 是否满足非共价条件 | `true` |
| `occurrences[].complete_heavy_atoms` | 可选布尔值或 null | 沉积坐标是否包含模板的全部重原子 | `false` |
| `occurrences[].template_compatible` | 可选布尔值或 null | 模板是否满足 Stage3 图结构和对称性条件 | `true` |
| `occurrences[].language_compatible` | 可选布尔值或 null | 身份是否与 Stage3 语言表征一致 | `true` |
| `assignments` | 长度 A 的列表 | A 为已判定归属的预测区域数，每项记录一个区域的归属 | 见上例，A = 1 |
| `assignments[].blob_id` | 非负整数 | Find 的预测区域编号 | `17` |
| `assignments[].occurrence_id` | 非负整数或 null | 所属拷贝的 `occurrences[].occurrence_id`，背景为 `null` | `0` |

归属未知的区域省略对应的 `assignments` 项。`narrow` 在所选 SMILES 与归属拷贝的身份相同时，按表中从 `eligible` 到 `language_compatible` 的顺序检查资格。第一个 `false` 决定跳过原因，缺项或 `null` 表示该条件未知，继续检查下一项。上例自动选 `CCO` 时，跳过原因为 `reference_failed_complete_heavy_atoms`。另一个乙醇拷贝可以有不同资格。

## Find 和 Match 的逐区域预测

### `find.json`

| 字段 | 类型或形状 | 编号、单位和实例 |
| --- | --- | --- |
| `sample_id`、`identity_smiles` | 字符串、字符串列表 (I,) | 如 `"demo-ca2"`、`["CCO","CCN"]`，在同一样本的四个阶段结果中共用。 |
| `seed` | 整数 | 本次 Find 种子，如 7。 |
| `candidates` | 对象列表，长度 B | 每项是一个原区域。零区域保存 `[]`。 |
| `blob_id` | 非负整数 | 如 17，连接完整来源区域，样本内唯一。 |
| `centered_box_index` | 非负整数 | 如 4，在局部特征 NPZ 的同名数组中查找编号 4 的位置。 |
| `center_xyz_A` | 数值列表 (3,) | 如 `[12.5,-3.0,8.25]`，整个来源区域所有体素中心在实验世界坐标中的 XYZ 均值，单位 Å。 |
| `selected` | 布尔 | 如 `true`，Find 的 Gaussian 选择标志；两组都进入考察。 |
| `probability_mean` | 数值，[0,1] | 如 0.74，完整来源区域的原平均概率。 |
| `identity_scores`、`foreground_probability` | Find 阶段均为 null | 如 `"identity_scores":null`，此时还没有 Match 预测。 |

### `match.json`

交接已有预测时，调用者须选取同一样本。每个区域的旧得分按其自身文件中的 `identity_smiles` 解释，再与本次样本的规范身份对齐。例如旧轴 `["CCN","CCO"]` 的 `[0.7,1.2]` 转到 `["CCO","CCN"]` 时为 `[1.2,0.7]`。旧身份归并到同一身份时取其最高有效分数，全部缺失则为 `null`。

Find 已执行且 B=0 时，整轮推理可写出空的 Find、Match 和计划候选列表，以及 `build.json` 的 `results=[]`。阶段错误会使进程返回非零状态，已完成阶段的文件仍可查看。

沿用上述顶层和区域字段，`seed` 改为 Match 种子，如 11。以下三项描述身份和前景。对空身份轴可以有 `identity_scores=[]`、`predicted_smiles=null`。

| 字段 | 类型或形状 | 含义与示例 |
| --- | --- | --- |
| `identity_scores` | 长度 I 的列表，每项为有限数值或 null | 如 `[1.2,0.7]`，按原身份轴 `["CCO","CCN"]` 顺序给出的 logits。缺少向量或无法构图的身份保留原位，得分为 `null`。 |
| `predicted_smiles` | 字符串或 null | 取第一个得分最高的候选，如 `"CCO"`。 |
| `foreground_probability` | 数值，[0,1] 或 null | 如 0.31，Match 对配体前景的概率，独立于身份得分。 |

## 对接计划和实际提交

### `decisions.json`

顶层是长度 M 的列表，M 为有人工决定的区域数，每个区域一项。构造示例为 `[{"blob_id":17,"action":"dock","smiles":"CCN","engine":"dg","seed":123}]`。

| 每项字段 | 类型 | 含义 | 构造示例 |
| --- | --- | --- | --- |
| `blob_id` | 非负整数 | Find 的原区域编号 | `17` |
| `action` | 字符串 | `skip` 跳过，`dock` 指定对接，`auto` 沿用自动选择 | `"dock"` |
| `smiles` | 可选字符串 | `dock` 手选的精确身份，省略时取 Match 预测。仅有 Find 结果时可直接指定。新身份从 SMILES 构图，原 Match 得分轴保持原顺序 | `"CCN"` |
| `engine` | 当前固定字符串 | 首版的计划固定写入 DG-PocketXMol。决定文件可省略此字段 | `"dg"` |
| `seed` | 可选整数 | `auto` 或 `dock` 的单请求采样种子，范围为 0 到 2³²−1 | `123` |

### `plan.json`

顶层有共同的 `sample_id` 和 `identity_smiles`。`scope` 是范围（如 `"wide"`），`k` 是样本内非负预算（如 1，`null` 不限），`seed` 是 Build 基础种子（如 0），`foreground_threshold` 是前景概率下限（如 0.4，`null` 关闭）。`candidates` 保留全部 B 个区域，按考察顺序排列。以下是与 README 输入配套的完整构造候选。

```json
{"blob_id":17,"centered_box_index":4,"center_xyz_A":[12.5,-3.0,8.25],"selected":true,"probability_mean":0.74,"identity_scores":[1.2,0.7],"foreground_probability":0.31,"predicted_smiles":"CCO","probability_rank":1,"inspection_rank":1,"chosen_smiles":"CCN","manual_input_smiles":"CCN","selection_source":"manual","plan_action":"dock","reason":"manual_dock","engine":"dg","sampling_seed":123,"planned_docking_rank":1}
```

| 增加的区域字段 | 类型或形状 | 解释与边界例 |
| --- | --- | --- |
| `probability_rank` | 正整数 | 原平均概率降序名次，从 1 起，同分取较小 `blob_id`。 |
| `inspection_rank` | 正整数 | 先考察 `selected=true` 的一组，再考察其余区域，组内沿用概率名次。原概率第 3 可是考察第 1。 |
| `manual_input_smiles` | 字符串或 null | 用户手选的原文，如 `"C[C@H](F)Cl"`；未指定时为 `null`。实际模型输入看 `chosen_smiles` |
| `chosen_smiles` | 字符串或 null | 实际选给 Build 的身份，如手选 `"CCN"`。手动跳过或无身份时为空，预算用完仍保留所选身份。 |
| `selection_source` | 字符串 | 如 `"manual"`，自动为 `"automatic"`。仅用 `auto` 覆盖种子仍属自动。 |
| `plan_action` | 字符串 | `"dock"` 安排一次请求，`"skip"` 保留该区域但不安排。 |
| `reason` | 字符串 | 手动为 `manual_skip` 或 `manual_dock`。自动通过为 `eligible`，自动跳过为 `no_identity`、`below_foreground_threshold`、`known_non_small_molecule` 或 `reference_failed_` 加第一个失败的资格名。预算耗尽为 `budget_exhausted`。 |
| `engine`、`sampling_seed` | 字符串、整数，或均为 null | 如 `"dg"` 和 123。缺少单区域种子时为 `(seed + blob_id) % 2**32`。 |
| `planned_docking_rank` | 正整数或 null | 从 1 起，仅对安排的请求连续编号。跳过者为 `null`，引擎和种子也为 `null`。 |
| `reference` | 可选实例字典或 null | 仅附本次窄范围读取的归属，字段与 `reference.json` 的 `occurrences[]` 相同，例如 `occurrence_id=0`、`smiles="CCO"`。给 Build 的执行参数只含编号、中心、所选身份、引擎、种子和计划号。 |

自动选择依次看是否有身份、可选前景阈值、类别和窄范围资格。已知前景概率低于阈值才筛除。`wide` 和 `narrow` 接受 `small_molecule` 与未知类别，`all` 允许全部类别。人工指定对接优先于这些自动条件，仍按相同考察次序消耗预算。

### `build.json`

顶层为共同的 `sample_id`、`identity_smiles` 和 `results`。`results` 是长度 R 的列表，R 为实际提交的对接请求数，包括已提交后失败的请求。以下为一个请求失败的完整构造记录。

```json
{"blob_id":17,"docking_rank":1,"status":"failed","pose_file":null,"top1_sample_index":null,"failure_reason":"no valid pose"}
```

| `results[]` 子字段 | 类型 | 含义 | 构造示例 |
| --- | --- | --- | --- |
| `blob_id` | 非负整数 | 本次请求对应的计划区域编号 | `17` |
| `docking_rank` | 正整数 | 实际提交次序，从 1 起 | `1` |
| `status` | 字符串 | `success` 表示姿态全部成功，`partial` 表示部分成功，`failed` 表示全部失败 | `"failed"` |
| `pose_file` | 路径字符串或 null | 成功姿态的 SDF 路径，相对 `build/`，无成功姿态时为 `null` | `"17/poses.sdf"` |
| `top1_sample_index` | 非负整数或 null | 得分最高的姿态在原采样中的编号，无可用评分时为 `null` | `9` |
| `failure_reason` | 字符串或 null | 请求的失败原因，成功时为 `null` | `"no valid pose"` |

#### `build/17/poses.sdf`

每个成功姿态保留原实验世界 XYZ 坐标，单位 Å，SDF 分子位置从 0 起。

#### `build/17/candidates.json`

顶层是长度 P 的列表，P 为该请求的采样姿态数。每项的 `sample_index` 是从 0 起的采样编号，`sdf_index` 是成功姿态在 `poses.sdf` 中从 0 起的位置，失败时填 `null`。例如原采样 9 可位于成功 SDF 的第 1 个零基位置，失败姿态没有 SDF 位置。

#### `build/17/ranking.json`

`candidate_metrics` 是长度 P 的列表，每项对应一个采样姿态。`candidate_metrics[].sample_index` 是原采样编号（如 `9`），`candidate_metrics[].self_ranking` 是综合评分（如 `2.4`），无法评分时为 `null`。评分为轨迹置信度 `cfd_traj + int(no_clashes) + int(stereo)`，后两项分别表示通过碰撞和立体检查。文件中的记录保留采样顺序及失败项。取前 k 个时，先选择具有综合评分的姿态，按分数降序排列，同分取原采样编号较小者，再取前 k 项。`top1_sample_index` 由此获得。

#### `build/17/confidence.npz`

S 为写入 SDF 的姿态数，N 为配体重原子数，T 为每个姿态的采样步数，H = N(N−1)/2 为不重复原子对数。各数组的第一个轴与 SDF 姿态同序。

| 数组 | 类型与形状 | 含义 |
| --- | --- | --- |
| `sample_index` | int64，(S,) | 每个成功姿态的原采样编号，如 `[1,9]` |
| `confidence_pos_traj` | float32，(S,N,T) | 每个姿态中每个原子在各步的位置置信度 |
| `confidence_pos` | float32，(S,N,1) | 最后一步每个原子的位置置信度 |
| `confidence_node` | float32，(S,N,1) | 最后一步每个原子的类别置信度 |
| `confidence_halfedge` | float32，(S,H,1) | 最后一步每个原子对的类别置信度，原子对按完全图上三角排列 |

姿态数量（论文为每请求 50 个，每姿态 100 步）与样本请求预算 K 独立。

## 运行记录

### `run.json`

`identity_inputs` 按 Infer 读到的 `sample.json.identities` 顺序保存 `input_smiles` 与规范 `smiles`。例如旧样本中的 `OCC` 可记为 `{"input_smiles":"OCC","smiles":"CCO"}`，无法解析时 `smiles` 为 `null`。从原始数据执行 `run` 时，此处的输入已由 Prepare 归并。用户最初提交的每种写法及其身份下标保存在 `prepared/prepare.json` 的 `smiles_inputs`。

`sample` 是与 [sample.json](#samplejson) 同结构的对象，其中地图、受体和身份表征路径均为绝对路径。`command` 是本次调用的参数对象。

| `command` 子字段 | 类型 | 含义 | 构造示例 |
| --- | --- | --- | --- |
| `sample`、`output` | 路径字符串 | 样本说明文件、本次输出目录的绝对路径 | `"/data/sample.json"`、`"/data/reviewed"` |
| `models`、`candidates`、`decisions`、`reference` | 路径字符串或 null | 依次指向模型目录、候选文件、人工决定文件、真值文件 | `"/data/match.json"` 是 `candidates` 的示例 |
| `command`、`until`、`scope` | 字符串 | 入口名称、结束阶段、自动选择范围 | `"infer"`、`"build"`、`"wide"` |
| `k`、`foreground_threshold` | 非负整数或 null、概率或 null | 实际请求上限、前景概率下限 | `1`、`0.4` |
| `find_seed`、`match_seed`、`build_seed` | 整数 | 三个阶段各自的随机种子 | `7`、`11`、`23` |

该文件在阶段调用前写出，可用于追溯未完成的运行。
