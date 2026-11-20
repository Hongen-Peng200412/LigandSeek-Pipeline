# 单样本产物与模型身份

`prepare_sample` 返回新样本的 `sample.json`。该文件引用的输入均位于 `--output` 下。以下两段目录树描述同一个物理目录 `prepared/`。编号 `<i>` 从 0 起，按第一个出现该规范身份的用户输入排序；只有成功得到向量的编号才有 NPZ。

```text
## <科学产物>
prepared/
├── sample.json                         # Infer 的单样本入口
├── receptor.cif                        # 本次使用的 mmCIF
├── exp.npy                             # 实验密度
├── exp.npz                             # 实验图的世界坐标网格
├── sim.npy                             # 受体在同一网格上的模拟密度
├── sim.npz                             # 模拟图网格、分辨率与工具版本
└── smi_ted_289m/<sample_id>/
    ├── candidate_<i>.npz               # 一个成功编码身份的 768 维向量
    └── results.jsonl                   # 每个规范身份的编码状态

## <其他文件>
prepared/
├── prepare.json                        # 原文、文件来源、执行模式和工具版本（运行诊断）
├── cryoatom2.log                       # 实际运行 CryoAtom2 时的输出（运行诊断）
└── smi_ted.log                         # 编码器命令、耗时和输出（运行诊断）
```

## 样本与候选身份

### sample.json

下面是两个输入 SMILES 对应同一个模型身份的完整构造示例。

```json
{
  "sample_id": "demo",
  "experimental_map": "exp.npy",
  "simulated_map": "sim.npy",
  "receptor_cif": "receptor.cif",
  "receptor_source": "ca2",
  "map_convention": "mrc",
  "identities": [
    {"smiles": "CCO", "embedding": "smi_ted_289m/demo/candidate_0.npz"}
  ]
}
```

`sample_id` 为小写样本标识。`receptor_source` 决定 Find 的受体校准，`deposited` 对应已知受体，`ca2` 对应预测受体。`map_convention` 记录原始实验图的坐标解释，两个 NPY 均已落在模型所用网格上。

`identities` 的列表下标是 Match 得分的身份下标。RDKit 用 canonical、non-isomeric SMILES 表示身份。例如 `C[C@H](O)C(=O)O` 和 `C[C@@H](O)C(=O)O` 共用 `CC(O)C(=O)O`，双键方向和同位素标记也不区分。可解析的身份即使编码失败仍占一个位置，例如 `{"smiles":"CCO","embedding":null}` 表示这个槽没有可用向量。完全无法解析的原文只出现在 `prepare.json` 的来源映射中。

### receptor.cif

这是一份世界 XYZ 坐标以 Å 计的 mmCIF。提供已有受体时复制用户给定的 CIF；重新建模时取 CryoAtom2 最终过滤后的 CIF。模拟图只使用该 CIF 首个模型中按原规则选定的 `ATOM` 重原子。Infer 独立按 CIF 的聚合物实体识别受体特征，CA2 未写实体表时按 `ATOM` 记录识别。

## 密度和坐标

### exp.npy 与 sim.npy

两个文件均为 `float32 (1,Z,Y,X)`，形状完全相同。第 0 轴是一个密度通道，后三轴以 Z、Y、X 排列。例如 `sim[0,2,3,4]` 是索引 XYZ 为 `[4,3,2]` 的受体模拟密度，`exp[0,2,3,4]` 是同一点的实验密度。实际 Z、Y、X 长度由原始样本决定。

### exp.npz

| 字段 | 数据类型、单位及含义 | 构造示例 |
| --- | --- | --- |
| `origin` | `float32 (3,)`，世界 XYZ 网格下角点，Å | `[-0.5,-0.5,-0.5]` |
| `voxel_size` | `float32 (3,)`，沿 XYZ 的真实间距，Å/体素 | `[1.0,1.0,1.0]` |
| `canonical_shape_zyx` | `int32 (3,)`，密度数组后三维的大小 | `[96,100,104]` |
| `schema_version` | `uint16` 标量，几何格式版本 | `2` |

体素中心按 `origin + ([x,y,z]+0.5)*voxel_size` 定位。在表示例中，索引 `[4,3,2]` 对应世界坐标 `[4,3,2]` Å。标准 MRC 的非零 ORIGIN 表示首个采样点；ORIGIN 为零时用按 MAPC/MAPR/MAPS 还原至 XYZ 的 NSTART 与体素间距计算首点。原图重采样后使用新的实际间距。历史 `adaligand` 路径沿用原始导入与 Chimera 写图的坐标约定。

### sim.npz

`origin`、`voxel_size`、`canonical_shape_zyx`、`schema_version` 与 `exp.npz` 对应字段相同。`resolution` 为 `float32` 标量，记录生成模拟图所用的 Å 数，例如 `3.2`。`chimera_version` 为字符串标量，保存可执行程序报告的实际版本。

## 分子向量

### smi_ted_289m/&lt;sample_id&gt;/candidate_&lt;i&gt;.npz

`embedding` 为 `float32 (768,)`，是一种规范身份的整体分子向量。标量 `candidate_id` 为 `int32`，即上述身份下标。其余标量是字符串。`pdb_id` 保存样本名，例如 `"demo"`，新样本无需拥有沉积 PDB 编号。`object_key` 是该候选的化学身份，在公开样本中例如 `"CCO"`。`model_name` 为 `"SMI-TED Light 289M"`。`prepared_smiles` 为送入官方编码器的规范 SMILES。`model_smiles` 为官方规范化函数得到的 tokenizer 输入，未能确认时以空字符串表示。

### smi_ted_289m/&lt;sample_id&gt;/results.jsonl

每个规范候选一行 JSON，按 `candidate_id` 递增。下面给出失败记录的构造示例。

```json
{"pdb_id":"demo","candidate_id":0,"type_tag":null,"status":"model_failed","prepared_smiles":"CCO","model_smiles":"CCO","output_file":null,"error":"encoder failed","diagnostics":{"official_normalization_success":true}}
```

`type_tag` 在盲编码时为 `null`，未知类别不影响向量计算。`status` 为 `encoded` 时，`output_file` 指本 JSONL 所在目录中的向量文件名，例如 `"candidate_0.npz"`，`error` 为 `null`。编码器失败时保存错误文本。没有输入身份的实例记录可为 `no_smiles`；公开 Prepare 已将无法解析的原文记在输入映射中。

`diagnostics` 来自原编码器包装。规范化诊断含调用字符串、规范化结果、成功值和错误。分词诊断含加入起止 token 前后的长度、模型上限、是否截断、截断范围内的前缀、未收录 token 列表和分词重组结果。例如 `unsupported_tokens: []` 表示本次拆出的全部 token 在词表内。成功写出向量时还有 `all_finite`、`nan_count`、`positive_inf_count`、`negative_inf_count`。例如 `all_finite: false` 时保留原向量供追溯，`sample.json` 将对应的 `embedding` 写为 `null`。

## 输入来源与执行记录

### prepare.json

| 字段 | 含义 | 构造示例 |
| --- | --- | --- |
| `input`、`experimental_map`、`models`、`smi_ted_model` | 准备时的输入说明、原图、模型资产根、SMI-TED 目录的路径；样本推理从 `sample.json` 读取科学输入。 | `"/work/demo/raw.json"` |
| `smiles_inputs` | 按用户原始输入顺序记录 `input_smiles` 与零基 `identity_index`。相同规范身份共享下标，解析失败时为 null。 | `[{"input_smiles":"OCC","identity_index":0},{"input_smiles":"CCO","identity_index":0},{"input_smiles":"bad?","identity_index":null}]` |
| `receptor` | 实际来源、模式、输入 CIF、FASTA 路径、实际参数与版本。命令内的临时路径只表明已用参数，正式受体在 `receptor.cif`。 | `{"source":"deposited","mode":"supplied","input_cif":"/work/r.cif","sequences":{},"command":null,"version":null}` |
| `density` | 受体模拟的实验分辨率、采用依据、网格约定、工具版本与用于模拟的原子数。 | `{"resolution":3.0,"resolution_source":"default","map_convention":"mrc","chimera":"chimera production version 1.19","atom_count":814}` |

`receptor.mode` 对 CA2 实际建模取 `with_sequences` 或 `without_sequences`。例如仅传 RNA 时，`sequences` 可为 `{"rna":"/work/rna.fa"}`。`density.resolution_source` 为 `user`、`emdb.final_reconstruction`、`rcsb.resolution_combined`、`mrc_label` 或 `default`。

### cryoatom2.log

只有本次调用 CryoAtom2 时才生成，记录命令、版本、标准输出和错误输出，末行记录耗时秒数和返回码。成功时 `prepare.json` 也记录 `receptor.elapsed_seconds`，如 `526.2`。

### smi_ted.log

保存实际编码命令、标准输出和错误输出，末行记录耗时秒数和返回码。例如 `{"elapsed_seconds":13.4,"returncode":0}`。工具或某个准备阶段失败时，仅保留相邻的 `<output>.unfinished-*` 目录，其中 `cryoatom2_work/` 或 `chimera_work/` 可保存尚未成功返回的程序原件，包括 Chimera 的 `molmap.stderr.log`。
